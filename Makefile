##########################
# Bootstrapping variables
##########################

# TODO: move to config file and parse with jq
export STAGE ?= "dev"
export STACK_NAME ?= "gfe-db-3"
export NEO4J_USERNAME ?= "neo4j"
export NEO4J_PASSWORD ?= "gfedb"

target:
	$(info ${HELP_MESSAGE})
	@exit 0

# init: ##=> Install OS deps and dev tools
# 	$(info [*] Bootstrapping CI system...)
# 	@$(MAKE) _install_os_packages

deploy: ##=> Deploy services
	$(info [*] Deploying...)
	bash scripts/deploy.sh ${STAGE} ${STACK_NAME} ${NEO4J_USERNAME} ${NEO4J_PASSWORD}

load: ##=> Load an IMGT/HLA release version
	$(info [*] Loading ${release})

#	# Get state machine arn

	aws stepfunctions start-execution \
		--state-machine-arn

# delete: ##=> Delete services
# 	aws cloudformation delete-stack --stack-name ${STACK_NAME}-${STAGE}

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

	These variables are automatically filled at CI time except STRIPE_SECRET_KEY
	If doing a dirty/individual/non-ci deployment locally you'd need them to be set

	STAGE: "dev"
		Description: Feature branch name used as part of stacks name
	STACK_NAME: "gfe-db"
		Description: Stack Name already deployed; used for dirty/individual deployment

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