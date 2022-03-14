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
export NEO4J_AMI_ID ?= ami-0e1324ddfc4d086bb # Requires subscription through AWS Marketplace
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


target:
	$(info ${HELP_MESSAGE})
	@exit 0

# TODO: Update email and name for Submitter node
deploy: logs.purge check-env ##=> Deploy services
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


check-env:
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

	ECR_BASE_URI: "$${ECR_BASE_URI}"
		Description: Base URI for AWS Elastic Container Registry

	BUILD_REPOSITORY: "$${BUILD_REPOSITORY}"
		Description: Name of the ECR repository for the build service

	LOAD_REPOSITORY: "$${LOAD_REPOSITORY}"
		Description: Name of the ECR repository for the load service

	Common usage:

	...::: Deploy all CloudFormation based services :::...
	$ make deploy

	...::: Run the StepFunctions State Machine to load Neo4j :::...
	$ make run release=3450 align=False kir=False mem_profile=False limit=1000

	...::: Delete all CloudFormation based services and data :::...
	$ make delete

endef
