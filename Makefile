##########################
# Bootstrapping variables
##########################

# TODO: move to config file and parse with jq, or use SSM parameter store
export STAGE ?= dev
export APP_NAME ?= gfe-db-3
export NEO4J_USERNAME ?= neo4j
export NEO4J_PASSWORD ?= gfedb
export REGION ?= us-east-1

target:
	$(info ${HELP_MESSAGE})
	@exit 0

# init: ##=> Install OS deps and dev tools
# 	$(info [*] Bootstrapping CI system...)
# 	@$(MAKE) _install_os_packages

deploy: ##=> Deploy resources
	$(info [*] Deploying...)
	@bash scripts/deploy.sh ${STAGE} ${APP_NAME} ${REGION} ${NEO4J_USERNAME} ${NEO4J_PASSWORD}

run: ##=> Load an IMGT/HLA release version; release=3450 align=False kir=False mem_profile=False limit=1000
	$(info [*] Starting StepFunctions execution for release $(release))

	# TODO: Add validation for positional arguments: release, align, kir, mem_profile, limit
	@aws stepfunctions start-execution \
	 	--state-machine-arn $$(aws ssm get-parameter \
	 		--name "/${APP_NAME}/${STAGE}/${REGION}/UpdatePipelineArn" | jq -r '.Parameter.Value') \
	 	--input "{\"params\":{\"environment\":{\"RELEASES\":\"$(release)\",\"ALIGN\":\"False\",\"KIR\":\"False\",\"MEM_PROFILE\":\"False\",\"LIMIT\":\"$(limit)\"}}}" > dev/null

delete: ##=> Delete resources
	$(info [*] Deleting resources...)
	@bash scripts/delete.sh ${STAGE} ${APP_NAME} ${REGION}

# export.parameter:
# 	$(info [+] Adding new parameter named "${NAME}")
# 	aws ssm put-parameter \
# 		--name "$${NAME}" \
# 		--type "String" \
# 		--value "$${VALUE}" \
# 		--overwrite

#############
#  Helpers  #
#############

# _install_os_packages:
# 	$(info [*] Installing jq...)
# 	yum install jq -y
# 	$(info [*] Upgrading Python SAM CLI and CloudFormation linter to the latest version...)
# 	python3 -m pip install --upgrade --user cfn-lint aws-sam-cli
# 	npm -g install aws-cdk

define HELP_MESSAGE

	Environment variables:

	STAGE: "dev"
		Description: Feature branch name used as part of stacks name
	APP_NAME: "gfe-db"
		Description: Stack Name already deployed

	Common usage:

	...::: Deploy all CloudFormation based services :::...
	$ make deploy

endef

# TODO: Makefile structure
# - deploy
# 	- setup
# 	- nested cfn (deploy.pipeline, deploy.database)
# 	- build, tag, push images to ECR (deploy.build, deploy.load)
# - load (version number): Run gfe-db (invoke Step Functions), example: make load 3450
# - delete
# 	- delete data in buckets