##########################
# Bootstrapping variables
##########################

# TODO: move to config file and parse with jq, or use SSM parameter store
export STAGE ?= dev
export APP_NAME ?= gfe-db
export NEO4J_USERNAME ?= neo4j
export NEO4J_PASSWORD ?= gfedb
export REGION ?= us-east-1

target:
	$(info ${HELP_MESSAGE})
	@exit 0


deploy: deploy.ecr
	$(info [*] Deploying resources...)

deploy.ecr: deploy.cfn
	$(info [*] Pushing Docker images to ECR...)
	@aws ecr get-login-password \
		--region ${REGION} | docker login \
			--username AWS \
			--password-stdin $$(aws sts get-caller-identity --query Account --output text).dkr.ecr.${REGION}.amazonaws.com

	@docker build -t ${STAGE}-${APP_NAME}-build-service build/
	@docker tag ${STAGE}-${APP_NAME}-build-service:latest $$(aws sts get-caller-identity --query Account --output text).dkr.ecr.${REGION}.amazonaws.com/${STAGE}-${APP_NAME}-build-service:latest
	@docker push $$(aws sts get-caller-identity --query Account --output text).dkr.ecr.${REGION}.amazonaws.com/${STAGE}-${APP_NAME}-build-service:latest

	@docker build -t ${STAGE}-${APP_NAME}-load-service load/
	@docker tag ${STAGE}-${APP_NAME}-load-service:latest $$(aws sts get-caller-identity --query Account --output text).dkr.ecr.${REGION}.amazonaws.com/${STAGE}-${APP_NAME}-load-service:latest
	@docker push $$(aws sts get-caller-identity --query Account --output text).dkr.ecr.${REGION}.amazonaws.com/${STAGE}-${APP_NAME}-load-service:latest

deploy.cfn:
	$(info [*] Deploying...)
	@bash scripts/deploy.sh ${STAGE} ${APP_NAME} ${REGION} ${NEO4J_USERNAME} ${NEO4J_PASSWORD}

run: ##=> Load an IMGT/HLA release version; release=3450 align=False kir=False mem_profile=False limit=1000
	$(info [*] Starting StepFunctions execution for release $(release))

#	@# TODO: Add validation for positional arguments: release, align, kir, mem_profile, limit
	@echo "Execution running:"
	@aws stepfunctions start-execution \
	 	--state-machine-arn $$(aws ssm get-parameter --name "/${APP_NAME}/${STAGE}/${REGION}/UpdatePipelineArn" | jq -r '.Parameter.Value') \
	 	--input "{\"params\":{\"environment\":{\"RELEASES\":\"$(release)\",\"ALIGN\":\"False\",\"KIR\":\"False\",\"MEM_PROFILE\":\"False\",\"LIMIT\":\"$(limit)\"}}}" | jq '.executionArn'

delete: ##=> Delete resources
	$(info [*] Deleting resources...)
	@bash scripts/delete.sh ${STAGE} ${APP_NAME} ${REGION}

define HELP_MESSAGE

	Environment variables:

	STAGE: "${STAGE}"
		Description: Feature branch name used as part of stacks name
	APP_NAME: "${APP_NAME}"
		Description: Stack Name already deployed
	REGION: "${REGION}":
		Description: AWS region for deployment
	NEO4J_USERNAME: "${NEO4J_USERNAME}"
		Description: Neo4j username
	NEO4J_PASSWORD: "${NEO4J_PASSWORD}"
		Description: Neo4j password

	Common usage:

	...::: Deploy all CloudFormation based services :::...
	$ make deploy

	...::: Run the StepFunctions State Machine to load Neo4j :::...
	$ make run release=3450 align=False kir=False mem_profile=False limit=1000

	...::: Delete all CloudFormation based services and data :::...
	$ make delete

endef
