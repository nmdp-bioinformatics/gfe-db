##########################
# Bootstrapping variables
##########################

# Application specific environment variables
include .env
export

# Base settings, these should almost never change
export AWS_ACCOUNT ?= $(shell aws sts get-caller-identity --query Account --output text)

# TODO: Application Configuration, can move to JSON
export ROOT_DIR ?= $(shell pwd)
export DATABASE_DIR ?= ${ROOT_DIR}/${APP_NAME}/database
export LOGS_DIR ?= $(shell echo "${ROOT_DIR}/logs")
export CFN_LOG_PATH ?= $(shell echo "${LOGS_DIR}/cfn/logs.txt")
export PURGE_LOGS ?= false

# TODO move these to a config file
export NEO4J_AMI_ID ?= ami-04aa5da301f99bf58 # Bitnami Neo4j, requires subscription through AWS Marketplace
export DATABASE_VOLUME_SIZE ?= 50
# TODO: Add TRIGGER_SCHEDULE variable
# TODO: Add BACKUP_SCHEDULE variable

# Resource identifiers
export DATA_BUCKET_NAME ?= ${STAGE}-${APP_NAME}-${AWS_ACCOUNT}-${AWS_REGION}
export ECR_BASE_URI ?= ${AWS_ACCOUNT}.dkr.ecr.${AWS_REGION}.amazonaws.com
export BUILD_REPOSITORY ?= ${STAGE}-${APP_NAME}-build-service
# export SUBDOMAIN ?= ${STAGE}-${APP_NAME}
export INSTANCE_ID ?= $(shell aws ssm get-parameters \
		--names "/${APP_NAME}/${STAGE}/${AWS_REGION}/Neo4jDatabaseInstanceId" \
		--output json \
		| jq -r '.Parameters | map(select(.Version == 1))[0].Value')

# S3 paths
export PIPELINE_STATE_PATH ?= config/IMGTHLA-repository-state.json
export PIPELINE_PARAMS_PATH ?= config/pipeline-input.json
export FUNCTIONS_PATH ?= ${APP_NAME}/pipeline/functions

# App state values
export INSTANCE_STATE ?= $(shell aws ec2 describe-instance-status | jq -r '.InstanceStatuses[] | select(.InstanceId | contains("i-0ea29a765388720a8")).InstanceState.Name')
export NEO4J_ENDPOINT ?= $(shell aws ssm get-parameters \
	--names "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/Neo4jDatabaseEndpoint" \
	| jq -r '.Parameters | map(select(.Version == 1))[0].Value')

target:
	$(info ${HELP_MESSAGE})
	@exit 0

# TODO: Update email and name for Submitter node
deploy: logs.purge check.env ##=> Deploy services
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying ${APP_NAME} to ${AWS_ACCOUNT}" 2>&1 | tee -a ${CFN_LOG_PATH}
	$(MAKE) infrastructure.deploy
	$(MAKE) database.deploy
	$(MAKE) pipeline.deploy
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
ifndef AWS_REGION
$(error AWS_REGION is not set. Please add AWS_REGION to the environment variables.)
endif
ifndef AWS_PROFILE
$(error AWS_PROFILE is not set. Please select an AWS profile to use.)
endif
ifndef GITHUB_PERSONAL_ACCESS_TOKEN
$(error GITHUB_PERSONAL_ACCESS_TOKEN is not set. Please add GITHUB_PERSONAL_ACCESS_TOKEN to the environment variables.)
endif
ifndef HOST_DOMAIN
$(error HOST_DOMAIN is not set. Please add HOST_DOMAIN to the environment variables.)
endif
ifndef ADMIN_EMAIL
$(error ADMIN_EMAIL is not set. Please add ADMIN_EMAIL to the environment variables.)
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
infrastructure.deploy:
	$(MAKE) -C ${APP_NAME}/infrastructure/ deploy

database.deploy:
	$(MAKE) -C ${APP_NAME}/database/ deploy

pipeline.deploy:
	$(MAKE) -C ${APP_NAME}/pipeline/ deploy

pipeline.jobs.deploy:
	$(MAKE) -C ${APP_NAME}/pipeline/jobs/ deploy

config.deploy:
	$(MAKE) -C ${APP_NAME}/pipeline/ config.deploy
	$(MAKE) -C ${APP_NAME}/database/ config.deploy

database.load:
	@echo "Confirm payload:" && \
	[ "$$align" ] && align="$$align" || align="False" && \
	[ "$$kir" ] && kir="$$kir" || kir="False" && \
	[ "$$limit" ] && limit="$$limit" || limit="" && \
	[ "$$releases" ] && releases="$$releases" || releases="" && \
	payload="{ \"align\": \"$$align\", \"kir\": \"$$kir\", \"limit\": \"$$limit\", \"releases\": \"$$releases\", \"mem_profile\": \"False\" }" && \
	echo "$$payload" | jq -r && \
	echo "$$payload" | jq > payload.json
	@echo "Run pipeline with this payload? [y/N] \c " && read ans && [ $${ans:-N} = y ]
	@function_name="${STAGE}"-"${APP_NAME}"-"$$(cat ${FUNCTIONS_PATH}/environment.json | jq -r '.Functions.InvokePipeline.FunctionConfiguration.FunctionName')" && \
	aws lambda invoke \
		--cli-binary-format raw-in-base64-out \
		--function-name "$$function_name" \
		--payload file://payload.json \
		response.json 2>&1

database.status:
	@echo "Current state: $$INSTANCE_STATE"
	 
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

# TODO account for http or https and whether or not EIP or DNS is being used
database.get-endpoint:
	@echo "http://$${NEO4J_ENDPOINT}:7473/browser/"

database.get-credentials:
	@secret_string=$$(aws secretsmanager get-secret-value --secret-id ${APP_NAME}-${STAGE}-Neo4jCredentials | jq -r '.SecretString') && \
	echo "Username: $$(echo $$secret_string | jq -r '.NEO4J_USERNAME')" && \
	echo "Password: $$(echo $$secret_string | jq -r '.NEO4J_PASSWORD')"

delete: # data=true/false ##=> Delete services
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting ${APP_NAME} in ${AWS_ACCOUNT}" 2>&1 | tee -a ${CFN_LOG_PATH}
	$(MAKE) pipeline.delete
	$(MAKE) database.delete
	$(MAKE) infrastructure.delete
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Finished deleting ${APP_NAME} in ${AWS_ACCOUNT}" 2>&1 | tee -a ${CFN_LOG_PATH}

# Delete specific stacks
infrastructure.delete:
	$(MAKE) -C ${APP_NAME}/infrastructure/ delete

database.delete:
	$(MAKE) -C ${APP_NAME}/database/ delete

pipeline.delete:
	$(MAKE) -C ${APP_NAME}/pipeline/ delete

pipeline.jobs.delete:
	$(MAKE) -C ${APP_NAME}/pipeline/jobs/ delete

# Administrative functions
get.data: #=> Download the build data locally
	@mkdir -p ${ROOT_DIR}/data
	@aws s3 cp --recursive s3://${DATA_BUCKET_NAME}/data/ ${ROOT_DIR}/data/

get.logs: #=> Download all logs locally
	@aws s3 cp --recursive s3://${DATA_BUCKET_NAME}/logs/ ${LOGS_DIR}/

# # TODO get pipeline execution status
# pipeline.status:

define HELP_MESSAGE

	Environment variables:

	STAGE: "${STAGE}"
		Description: Feature branch name used as part of stacks name

	APP_NAME: "${APP_NAME}"
		Description: Stack Name already deployed

	AWS_ACCOUNT: "${AWS_ACCOUNT}":
		Description: AWS account ID for deployment

	AWS_REGION: "${AWS_REGION}":
		Description: AWS region for deployment

	DATA_BUCKET_NAME "$${DATA_BUCKET_NAME}"
		Description: Name of the S3 bucket for data, config and logs

	Common usage:

	...::: Deploy all CloudFormation based services :::...
	$ make deploy

	...::: Deploy config files and scripts to S3 :::...
	$ make config.deploy

	...::: Run the StepFunctions State Machine to load Neo4j :::...
	$ make database.load releases=<version> align=<boolean> kir=<boolean> limit=<int>

	...::: Download CSV data from S3 :::...
	$ make get.data

	...::: Download logs from EC2 :::...
	$ make get.logs

	...::: Delete all CloudFormation based services and data :::...
	$ make delete

endef
