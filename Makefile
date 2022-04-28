##########################
# Bootstrapping variables
##########################

# Base settings, these should almost never change
export STAGE ?= dev
export APP_NAME ?= gfe-db
export AWS_ACCOUNT ?= $(shell aws sts get-caller-identity --query Account --output text)
export REGION ?= us-east-1

# TODO: Application Configuration, can move to JSON
export ROOT_DIR ?= $(shell pwd)
export LOGS_DIR ?= $(shell echo "${ROOT_DIR}/logs")
export CFN_LOG_PATH ?= $(shell echo "${LOGS_DIR}/cfn/logs.txt")
export PURGE_LOGS ?= false
# export NEO4J_AMI_ID ?= ami-0e1324ddfc4d086bb # Neo4j CE is discontinued; Requires subscription through AWS Marketplace
export NEO4J_AMI_ID ?= ami-04aa5da301f99bf58 # Bitnami Neo4j, requires subscription through AWS Marketplace
export DATABASE_VOLUME_SIZE ?= 50
# TODO: Add TRIGGER_SCHEDULE variable
# TODO: Add BACKUP_SCHEDULE variable

# Resource identifiers
export DATA_BUCKET_NAME ?= ${STAGE}-${APP_NAME}-${AWS_ACCOUNT}-${REGION}
export ECR_BASE_URI ?= ${AWS_ACCOUNT}.dkr.ecr.${REGION}.amazonaws.com
export BUILD_REPOSITORY ?= ${STAGE}-${APP_NAME}-build-service
export LOAD_REPOSITORY ?= ${STAGE}-${APP_NAME}-load-service
export PIPELINE_STATE_PATH ?= config/IMGTHLA-repository-state.json
export PIPELINE_PARAMS_PATH ?= config/pipeline-input.json
export FUNCTIONS_PATH ?= ${APP_NAME}/pipeline/functions

export INSTANCE_ID ?= $(shell aws ssm get-parameters \
		--names "/${APP_NAME}/${STAGE}/${REGION}/Neo4jDatabaseInstanceId" \
		--output json \
		| jq -r '.Parameters | map(select(.Version == 1))[0].Value')
export NEO4J_ENDPOINT=$(shell aws ssm get-parameters \
	--names "/$${APP_NAME}/$${STAGE}/$${REGION}/Neo4jDatabaseEndpoint" \
	| jq -r '.Parameters | map(select(.Version == 1))[0].Value')

# # Capture datetime of most recent parameter change (force refresh paramter references)
# export SSM_PARAM_MODIFIED ?= $(shell aws ssm describe-parameters \
# 	| jq -c '.Parameters[] | select(.Name | contains("/${APP_NAME}/${STAGE}/${REGION}/"))' \
# 	| jq -r '.LastModifiedDate' | sort -r | head -n 1)

target:
	$(info ${HELP_MESSAGE})
	@exit 0

# TODO: Update email and name for Submitter node
deploy: logs.purge check.env ##=> Deploy services
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying ${APP_NAME} to ${AWS_ACCOUNT}" 2>&1 | tee -a ${CFN_LOG_PATH}
	$(MAKE) deploy.infrastructure
	$(MAKE) deploy.database
	$(MAKE) deploy.pipeline
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Finished deploying ${APP_NAME}" 2>&1 | tee -a ${CFN_LOG_PATH}

logs.purge: logs.dirs
ifeq ($(PURGE_LOGS),true)
	@rm ${LOGS_DIR}/cfn/*.txt
endif

logs.dirs:
	@mkdir -p "${LOGS_DIR}/cfn" \
		"${LOGS_DIR}/pipeline/build" \
		"${LOGS_DIR}/pipeline/load" \
		"${LOGS_DIR}/database/bootstrap" || true


check.env: check.dependencies
ifndef AWS_PROFILE
$(error AWS_PROFILE is not set. Please select an AWS profile to use.)
endif
ifndef NEO4J_USERNAME
$(error NEO4J_USERNAME is not set.)
endif
ifndef NEO4J_PASSWORD
$(error NEO4J_PASSWORD is not set.)
endif
ifndef GITHUB_PERSONAL_ACCESS_TOKEN
$(error GITHUB_PERSONAL_ACCESS_TOKEN is not set.)
endif
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Found environment variables" 2>&1 | tee -a ${CFN_LOG_PATH}

check.dependencies:
	$(MAKE) check.dependencies.docker
	$(MAKE) check.dependencies.awscli
	$(MAKE) check.dependencies.samcli
	$(MAKE) check.dependencies.jq

check.dependencies.docker:
	@if ! docker info >/dev/null 2>&1; then \
		echo "**** Docker is not running. Please start Docker before deploying. ****" && \
		echo "**** Please refer to the documentation for a list of prerequisistes. ****" && \
		exit 1; \
	fi

check.dependencies.awscli:
	@if ! aws --version >/dev/null 2>&1; then \
		echo "**** AWS CLI not found. Please install AWS CLI before deploying. ****" && \
		echo "**** Please refer to the documentation for a list of prerequisistes. ****" && \
		exit 1; \
	fi

check.dependencies.samcli:
	@if ! sam --version >/dev/null 2>&1; then \
		echo "**** SAM CLI not found. Please install SAM CLI before deploying. ****" && \
		echo "**** Please refer to the documentation for a list of prerequisistes. ****" && \
		exit 1; \
	fi

check.dependencies.jq:
	@if ! jq --version >/dev/null 2>&1; then \
		echo "**** jq not found. Please install jq before deploying. ****" && \
		echo "**** Please refer to the documentation for a list of prerequisistes. ****" && \
		exit 1; \
	fi

# Deploy specific stacks
deploy.infrastructure:
	$(MAKE) -C gfe-db/infrastructure/ deploy

deploy.database:
	$(MAKE) -C gfe-db/database/ deploy

deploy.pipeline:
	$(MAKE) -C gfe-db/pipeline/ deploy

deploy.config:
	$(MAKE) -C gfe-db/pipeline/ deploy.config
	$(MAKE) -C gfe-db/database/ deploy.config

database.load:
	@echo "Confirm payload:" && \
	[ "$$align" ] && align="$$align" || align="False" && \
	[ "$$kir" ] && kir="$$kir" || kir="False" && \
	[ "$$limit" ] && limit="$$limit" || limit="" && \
	[ "$$releases" ] && releases="$$releases" || releases="" && \
	payload="{ \"align\": \"$$align\", \"kir\": \"$$kir\", \"limit\": \"$$limit\", \"releases\": \"$$releases\", \"mem_profile\": \"False\" }" && \
	echo "$$payload" | jq -r && \
	echo "$$payload" | jq > payload.json
	@echo -n "Run pipeline with this payload? [y/N] " && read ans && [ $${ans:-N} = y ]
	@function_name="${STAGE}"-"${APP_NAME}"-"$$(cat ${FUNCTIONS_PATH}/environment.json | jq -r '.Functions.InvokePipeline.FunctionConfiguration.FunctionName')" && \
	aws lambda invoke \
		--cli-binary-format raw-in-base64-out \
		--function-name "$$function_name" \
		--payload file://payload.json \
		response.json 2>&1

# TODO: fix
# database.status:
# 	@status=$$(aws ec2 describe-instance-status --instance-ids ${INSTANCE_ID} | jq -r '.InstanceStatuses[] | select(.InstanceId | contains("${INSTANCE_ID}")).CurrentState.Name') && \
# 	if [[ $$status = "" ]]; then status="stopped"; fi; \
# 	printf "Current state: $$status\n" && \
# 	[ "$$status" = "running" ] && printf "Neo4j is available at:\n$$NEO4J_ENDPOINT\n" || exit 0
	 

database.start:
	@echo "Starting $${APP_NAME} server..."
	@response=$$(aws ec2 start-instances --instance-ids ${INSTANCE_ID}) && \
	echo "Previous state: $$(echo "$$response" | jq -r '.StartingInstances[] | select(.InstanceId | contains("${INSTANCE_ID}")).PreviousState.Name')" && \
	echo "Current state: $$(echo "$$response" | jq -r '.StartingInstances[] | select(.InstanceId | contains("${INSTANCE_ID}")).CurrentState.Name')"

database.stop:
	@echo "Stopping $${APP_NAME} server..."
	@response=$$(aws ec2 stop-instances --instance-ids ${INSTANCE_ID}) && \
	echo "Previous state: $$(echo "$$response" | jq -r '.StoppingInstances[] | select(.InstanceId | contains("${INSTANCE_ID}")).PreviousState.Name')" && \
	echo "Current state: $$(echo "$$response" | jq -r '.StoppingInstances[] | select(.InstanceId | contains("${INSTANCE_ID}")).CurrentState.Name')"

delete: # data=true/false ##=> Delete services
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting ${APP_NAME} in ${AWS_ACCOUNT}" 2>&1 | tee -a ${CFN_LOG_PATH}
	$(MAKE) delete.pipeline
	$(MAKE) delete.database
	$(MAKE) delete.infrastructure
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Finished deleting ${APP_NAME} in ${AWS_ACCOUNT}" 2>&1 | tee -a ${CFN_LOG_PATH}

# Delete specific stacks
delete.infrastructure:
	$(MAKE) -C gfe-db/infrastructure/ delete

delete.database:
	$(MAKE) -C gfe-db/database/ delete

delete.pipeline:
	$(MAKE) -C gfe-db/pipeline/ delete

# Administrative functions
get.data: #=> Download the build data locally
	@mkdir -p ${ROOT_DIR}/data
	@aws s3 cp --recursive s3://${DATA_BUCKET_NAME}/data/ ${ROOT_DIR}/data/

get.logs: #=> Download all logs locally
	@aws s3 cp --recursive s3://${DATA_BUCKET_NAME}/logs/ ${LOGS_DIR}/

get.neo4j:
	@echo "http://$${NEO4J_ENDPOINT}:7474/browser/"

# TODO: finished administrative targets
# get.config:
# ifndef dir=""
# 	@aws s3 cp --recursive s3://${DATA_BUCKET_NAME}/data/ ${ROOT_DIR}/data/ 
# endif
# 	# @aws s3 cp --recursive s3://${DATA_BUCKET_NAME}/data/ $(dir)

# show.config:
# get.state:
# show.state:
# show.endpoint:

# run: ##=> Load an IMGT/HLA release version; make run release=3450 align=False kir=False mem_profile=False limit=1000
# 	$(info [*] Starting StepFunctions execution for release $(release))

# #	@# TODO: Add validation for positional arguments: release, align, kir, mem_profile, limit
# 	@echo "Execution running:"
# 	@aws stepfunctions start-execution \
# 	 	--state-machine-arn $$(aws ssm get-parameter --name "/${APP_NAME}/${STAGE}/${REGION}/UpdatePipelineArn" | jq -r '.Parameter.Value') \
# 	 	--input "{\"params\":{\"environment\":{\"RELEASES\":\"$(release)\",\"ALIGN\":\"False\",\"KIR\":\"False\",\"MEM_PROFILE\":\"False\",\"LIMIT\":\"$(limit)\"}}}" | jq '.executionArn'

define HELP_MESSAGE

	Environment variables:

	STAGE: "${STAGE}"
		Description: Feature branch name used as part of stacks name

	APP_NAME: "${APP_NAME}"
		Description: Stack Name already deployed

	AWS_ACCOUNT: "${AWS_ACCOUNT}":
		Description: AWS account ID for deployment

	REGION: "${REGION}":
		Description: AWS region for deployment

	DATA_BUCKET_NAME "$${DATA_BUCKET_NAME}"
		Description: Name of the S3 bucket for data, config and logs

	Common usage:

	...::: Deploy all CloudFormation based services :::...
	$ make deploy

	...::: Deploy config files and scripts to S3 :::...
	$ make deploy.config

	...::: Run the StepFunctions State Machine to load Neo4j :::...
	$ make database.load releases=<version> align=<boolean> kir=<boolean> limit=<int>

	...::: Download CSV data from S3 :::...
	$ make get.data

	...::: Download logs from EC2 :::...
	$ make get.logs

	...::: Display the Neo4j Browser endpoint URL :::...
	$ make get.neo4j

	...::: Delete all CloudFormation based services and data :::...
	$ make delete

endef
