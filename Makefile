##########################
# Bootstrapping variables
##########################
export STAGE ?= dev
export APP_NAME ?= gfe-db
export AWS_ACCOUNT ?= $(shell aws sts get-caller-identity --query Account --output text)
export REGION ?= us-east-1

export DATA_BUCKET_NAME ?= ${STAGE}-${APP_NAME}-${AWS_ACCOUNT}-${REGION}
export ECR_BASE_URI ?= ${AWS_ACCOUNT}.dkr.ecr.${REGION}.amazonaws.com
export BUILD_REPOSITORY ?= ${STAGE}-${APP_NAME}-build-service
export LOAD_REPOSITORY ?= ${STAGE}-${APP_NAME}-load-service

target:
	$(info ${HELP_MESSAGE})
	@exit 0

deploy: ##=> Deploy services
	$(info [*] Deploying all services to ${AWS_ACCOUNT}...)
	$(MAKE) deploy.infrastructure
	$(MAKE) deploy.database
	$(MAKE) deploy.pipeline
	@echo "Finished deploying ${APP_NAME}."

# Deploy specific stacks
deploy.infrastructure:
	$(MAKE) -C gfe-db/infrastructure/ deploy

deploy.database:
	$(MAKE) -C gfe-db/database/ deploy

deploy.pipeline:
	$(MAKE) -C gfe-db/pipeline/ deploy

delete: ##=> Delete services
	$(MAKE) delete.pipeline
	$(MAKE) delete.database
	$(MAKE) delete.infrastructure

# Delete specific stacks
delete.infrastructure:
	$(MAKE) -C gfe-db/infrastructure/ delete

delete.database:
	$(MAKE) -C gfe-db/database/ delete

delete.pipeline:
	$(MAKE) -C gfe-db/pipeline/ delete




deploy.ecr:
	$(info [*] Logging into ECR...)
	@aws ecr get-login-password \
		--region ${REGION} | docker login \
			--username AWS \
			--password-stdin $$(aws sts get-caller-identity --query Account --output text).dkr.ecr.${REGION}.amazonaws.com

#	# @$(info [*] Pushing build service image to ECR...)
	@docker build -t ${STAGE}-${APP_NAME}-build-service build/ && \
	docker tag ${STAGE}-${APP_NAME}-build-service:latest $$(aws sts get-caller-identity --query Account --output text).dkr.ecr.${REGION}.amazonaws.com/${STAGE}-${APP_NAME}-build-service:latest && \
	docker push $$(aws sts get-caller-identity --query Account --output text).dkr.ecr.${REGION}.amazonaws.com/${STAGE}-${APP_NAME}-build-service:latest

#	# @$(info [*] Pushing load service image to ECR...)
	@docker build -t ${STAGE}-${APP_NAME}-load-service load/ && \
	docker tag ${STAGE}-${APP_NAME}-load-service:latest $$(aws sts get-caller-identity --query Account --output text).dkr.ecr.${REGION}.amazonaws.com/${STAGE}-${APP_NAME}-load-service:latest && \
	docker push $$(aws sts get-caller-identity --query Account --output text).dkr.ecr.${REGION}.amazonaws.com/${STAGE}-${APP_NAME}-load-service:latest

deploy.cfn:
	$(info [*] Deploying...)
	@bash scripts/deploy.sh ${STAGE} ${APP_NAME} ${REGION}

run: ##=> Load an IMGT/HLA release version; make run release=3450 align=False kir=False mem_profile=False limit=1000
	$(info [*] Starting StepFunctions execution for release $(release))

#	@# TODO: Add validation for positional arguments: release, align, kir, mem_profile, limit
	@echo "Execution running:"
	@aws stepfunctions start-execution \
	 	--state-machine-arn $$(aws ssm get-parameter --name "/${APP_NAME}/${STAGE}/${REGION}/UpdatePipelineArn" | jq -r '.Parameter.Value') \
	 	--input "{\"params\":{\"environment\":{\"RELEASES\":\"$(release)\",\"ALIGN\":\"False\",\"KIR\":\"False\",\"MEM_PROFILE\":\"False\",\"LIMIT\":\"$(limit)\"}}}" | jq '.executionArn'

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
