##########################
# Bootstrapping variables
##########################

# Environment variables
# include .env # Optional, include STAGE and AWS_PROFILE
include .env.${STAGE}
export

SPLASH_FONT := slant

export AWS_ACCOUNT ?= $(shell aws sts get-caller-identity \
	--query Account \
	--output text)
export HAS_STAGE := $(shell aws ssm get-parameters \
		--names "/${APP_NAME}/${STAGE}/${AWS_REGION}/Stage" \
		--output json \
		| jq -r '.Parameters[0].Value')

export ROOT_DIR := $(shell pwd)
export DATABASE_DIR := ${ROOT_DIR}/${APP_NAME}/database
export INFRA_DIR := ${ROOT_DIR}/${APP_NAME}/infrastructure
export LOGS_DIR := $(shell echo "${ROOT_DIR}/logs")
export CFN_LOG_PATH := $(shell echo "${LOGS_DIR}/cfn/logs.txt")
export PURGE_LOGS := false

# conditionally required variable defaults
CREATE_VPC ?= false
USE_PRIVATE_SUBNET ?= false
CREATE_SSM_VPC_ENDPOINT ?=
CREATE_SECRETSMANAGER_VPC_ENDPOINT ?=
DEPLOY_NAT_GATEWAY ?=
DEPLOY_BASTION_SERVER ?=
VPC_ID ?=
PUBLIC_SUBNET_ID ?=
HOST_DOMAIN ?=
SUBDOMAIN ?=
PRIVATE_SUBNET_ID ?=
ADMIN_IP ?=
export DATABASE_VOLUME_SIZE ?= 64

# TODO move these to a config file
# TODO: Add TRIGGER_SCHEDULE variable
# TODO: Add BACKUP_SCHEDULE variable

# Resource identifiers
export DATA_BUCKET_NAME ?= ${STAGE}-${APP_NAME}-${AWS_ACCOUNT}-${AWS_REGION}
export ECR_BASE_URI := ${AWS_ACCOUNT}.dkr.ecr.${AWS_REGION}.amazonaws.com
export BUILD_REPOSITORY ?= ${STAGE}-${APP_NAME}-build-service
export EC2_KEY_PAIR_NAME := $${STAGE}-$${APP_NAME}-$${AWS_REGION}-neo4j-key
export INSTANCE_ID := $(shell aws ssm get-parameters \
	--names "/${APP_NAME}/${STAGE}/${AWS_REGION}/Neo4jDatabaseInstanceId" \
	--output json \
	| jq -r '.Parameters[0].Value')

# S3 paths
export PIPELINE_STATE_PATH := config/IMGTHLA-repository-state.json
export PIPELINE_PARAMS_PATH := config/pipeline-input.json
export FUNCTIONS_PATH := ${APP_NAME}/pipeline/functions

# TODO validate data types
# Required environment variables
REQUIRED_VARS := STAGE APP_NAME AWS_ACCOUNT AWS_REGION AWS_PROFILE SUBSCRIBE_EMAILS \
	GITHUB_REPOSITORY_OWNER GITHUB_REPOSITORY_NAME GITHUB_PERSONAL_ACCESS_TOKEN \
	ADMIN_EMAIL NEO4J_AMI_ID APOC_VERSION GDS_VERSION

BOOLEAN_VARS := CREATE_VPC USE_PRIVATE_SUBNET CREATE_SSM_VPC_ENDPOINT CREATE_SECRETSMANAGER_VPC_ENDPOINT \
	DEPLOY_NAT_GATEWAY DEPLOY_BASTION_SERVER

# stdout colors
# blue: runtime message, no action required
# green: parameter value message, no action required
# yellow: message to user, action required
# red: error message, action required
define blue
	@tput setaf 4
	@echo $1
	@tput sgr0
endef

define green
	@tput setaf 2
	@echo $1
	@tput sgr0
endef

define yellow
	@tput setaf 3
	@echo $1
	@tput sgr0
endef

define red
	@tput setaf 1
	@echo $1
	@tput sgr0
endef

target:
	$(info ${HELP_MESSAGE})
	@exit 0

test:
	$(MAKE) -C ${APP_NAME}/infrastructure/access-services/bastion-server/ service.authorize-security-group-ingress

splash-screen:
ifeq ($(SPLASH_FONT),slant)
	@echo "\033[0;34m                                            "
	@echo "\033[0;34m           ____                      __ __  "
	@echo "\033[0;34m   ____ _ / __/___              ____/ // /_ "
	@echo "\033[0;32m  / __ \`// /_ / _ \   ______   / __  // __ \\"
	@echo "\033[0;32m / /_/ // __//  __/  /_____/  / /_/ // /_/ /"
	@echo "\033[0;34m \__, //_/   \___/            \____//_____/ "
	@echo "\033[0;34m/____/                                      "
	@echo "\033[0;34m                                            "
	@echo "\033[0;34mCopyright © 2002-2023 National Marrow Donor Program. All rights reserved."
	@echo "\033[0;34m                                            \033[0m"
endif

env.print:
	@echo "\033[0;33mReview the contents of the .env file:\033[0m"
	@echo "+---------------------------------------------------------------------------------+"
	@awk '{ if (substr($$0, 1, 1) != "#") { line = substr($$0, 1, 76); if (length($$0) > 76) line = line "..."; printf "| %-79s |\n", line }}' .env.${STAGE}
	@echo "+---------------------------------------------------------------------------------+"
	@echo "\033[0;33mPlease confirm the above values are correct.\033[0m"

deploy: splash-screen logs.purge env.validate ##=> Deploy all services
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying ${APP_NAME} to ${AWS_ACCOUNT}" 2>&1 | tee -a ${CFN_LOG_PATH}
	$(MAKE) env.print
	@echo "Deploy stack to the \`${STAGE}\` environment? [y/N] \c " && read ans && [ $${ans:-N} = y ]
	$(MAKE) infrastructure.deploy 
	$(MAKE) database.deploy
	$(MAKE) pipeline.deploy
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Finished deploying ${APP_NAME}" 2>&1 | tee -a ${CFN_LOG_PATH}
	$(MAKE) options-screen

logs.purge: logs.dirs
ifeq ($(PURGE_LOGS),true)
	@rm ${LOGS_DIR}/cfn/*.txt
endif

logs.dirs:
	@mkdir -p "${LOGS_DIR}/cfn" \
		"${LOGS_DIR}/pipeline/build" \
		"${LOGS_DIR}/pipeline/load" \
		"${LOGS_DIR}/database/bootstrap" || true

check.dependencies:
	$(MAKE) check.dependencies.docker
	$(MAKE) check.dependencies.awscli
	$(MAKE) check.dependencies.samcli
	$(MAKE) check.dependencies.jq
	$(MAKE) check.dependencies.coreutils

check.dependencies.docker:
	@if docker info 2>&1 | grep -q 'Is the docker daemon running?'; then \
		echo "\033[0;31m**** Docker is not running. Please start Docker before deploying. ****\033[0m"; \
		echo "\033[0;31m**** Please refer to the documentation for a list of prerequisites. ****\033[0m"; \
		exit 1; \
	fi

check.dependencies.awscli:
	@if ! aws --version >/dev/null 2>&1; then \
		echo "\033[0;31m**** AWS CLI not found. Please install AWS CLI before deploying. ****\033[0m"; \
		echo "\033[0;31m**** Please refer to the documentation for a list of prerequisites. ****\033[0m"; \
		exit 1; \
	fi

check.dependencies.samcli:
	@if ! sam --version >/dev/null 2>&1; then \
		echo "\033[0;31m**** SAM CLI not found. Please install SAM CLI before deploying. ****\033[0m"; \
		echo "\033[0;31m**** Please refer to the documentation for a list of prerequisites. ****\033[0m"; \
		exit 1; \
	fi

check.dependencies.jq:
	@if ! jq --version >/dev/null 2>&1; then \
		echo "\033[0;31m**** jq not found. Please install jq before deploying. ****\033[0m"; \
		echo "\033[0;31m**** Please refer to the documentation for a list of prerequisites. ****\033[0m"; \
		exit 1; \
	fi

check.dependencies.coreutils:
	@if ! gdate --version >/dev/null 2>&1; then \
		echo "\033[0;31m**** GNU coreutils not found. Please install GNU coreutils before deploying. ****\033[0m"; \
		echo "\033[0;31m**** Please refer to the documentation for a list of prerequisites. ****\033[0m"; \
		exit 1; \
	fi

env.validate.stage:
	@res=$$(aws ssm get-parameters \
		--names "/${APP_NAME}/${STAGE}/${AWS_REGION}/Stage" \
		--output json \
		| jq -r '.Parameters[0].Value') && \
	[[ $$res = "null" ]] && echo "No deployed stage found" || echo "Found deployed stage: $$res" && \
	if [ "$${res}" = "null" ]; then \
		echo "\033[0;32m**** Starting new deployment. ****\033[0m"; \
	elif [ "$${res}" = "${STAGE}" ]; then \
		echo "\033[0;32m**** Found existing deployment for \`${STAGE}\` ****\033[0m"; \
	else \
		echo "\033[0;31m**** STAGE mismatch or bad credential configuration. ****\033[0m" && \
		echo "\033[0;31m**** Please refer to the documentation for a list of prerequisites. ****\033[0m" && \
		exit 1; \
	fi

env.validate.subdomain:
	@res=$$(aws route53 list-resource-record-sets --hosted-zone-id ${HOSTED_ZONE_ID} | \
		jq -r --arg fqdn "$$fqdn" '.ResourceRecordSets[] | select(.Name == $$fqdn) | .Name') && \
	[[ $$res = "" ]] && echo "\033[0;31mERROR: No Route53 domain found for $$fqdn\033[0m" && exit 1 || true

# env.validate.use-existing-vpc.vars:
# ifeq ($(VPC_ID),)
# 	$(call red, "\`VPC_ID\` must be set as an environment variable when \`CREATE_VPC\` is \`false\`")
# 	@exit 1
# else
# 	$(call green, "Found VPC_ID: ${VPC_ID}")
# endif
# ifeq ($(PUBLIC_SUBNET_ID),)
# 	$(call red, "\`PUBLIC_SUBNET_ID\` must be set as an environment variable when \`CREATE_VPC\` is \`false\`")
# 	@exit 1
# else
# 	$(call green, "Found PUBLIC_SUBNET_ID: ${PUBLIC_SUBNET_ID}")
# endif
# 	$(MAKE) env.validate.use-private-subnet.vars

env.validate.use-private-subnet.vars:
ifeq ($(USE_PRIVATE_SUBNET),true)
ifeq ($(CREATE_VPC),false)
ifeq ($(PUBLIC_SUBNET_ID),)
	$(call red, "\`PUBLIC_SUBNET_ID\` must be set as an environment variable when \`USE_PRIVATE_SUBNET\` is \`true\`")
	@exit 1
else
	$(call green, "Found PRIVATE_SUBNET_ID: ${PRIVATE_SUBNET_ID}")
endif
ifeq ($(PRIVATE_SUBNET_ID),)
	$(call red, "\`PRIVATE_SUBNET_ID\` must be set as an environment variable when \`USE_PRIVATE_SUBNET\` is \`true\`")
	@exit 1
else
	$(call green, "Found PRIVATE_SUBNET_ID: ${PRIVATE_SUBNET_ID}")
endif
endif
ifeq ($(ADMIN_IP),)
	$(call red, "\`ADMIN_IP\` must be set as an environment variable when \`USE_PRIVATE_SUBNET\` is \`true\`")
	@exit 1
endif
ifeq ($(DEPLOY_NAT_GATEWAY),)
	$(call red, "\`DEPLOY_NAT_GATEWAY\` must be set when \`USE_PRIVATE_SUBNET\` is \`true\`")
	@exit 1
endif
ifeq ($(DEPLOY_BASTION_SERVER),)
	$(call red, "\`DEPLOY_BASTION_SERVER\` must be set when \`USE_PRIVATE_SUBNET\` is \`true\`")
	@exit 1
endif
else ifeq ($(USE_PRIVATE_SUBNET),false)
	$(call blue, "**** This deployment uses a public subnet for Neo4j ****")
ifneq ($(DEPLOY_NAT_GATEWAY),)
	$(call red, "\`DEPLOY_NAT_GATEWAY\` must not be set when \`USE_PRIVATE_SUBNET\` is \`false\`")
	@exit 1
endif
ifneq ($(DEPLOY_BASTION_SERVER),)
	$(call red, "\`DEPLOY_BASTION_SERVER\` must not be set when \`USE_PRIVATE_SUBNET\` is \`false\`")
	@exit 1
endif
ifeq ($(HOST_DOMAIN),)
	$(call red, "\`HOST_DOMAIN\` must be set as an environment variable when \`USE_PRIVATE_SUBNET\` is \`false\`")
	@exit 1
else
	$(call green, "Found HOST_DOMAIN: ${HOST_DOMAIN}")
endif
ifeq ($(SUBDOMAIN),)
	$(call red, "\`SUBDOMAIN\` must be set as an environment variable when \`USE_PRIVATE_SUBNET\` is \`false\`")
	@exit 1
else
	$(call green, "Found SUBDOMAIN: ${SUBDOMAIN}")
endif
ifeq ($(HOSTED_ZONE_ID),)
	$(call red, "\`HOSTED_ZONE_ID\` must be set as an environment variable when \`USE_PRIVATE_SUBNET\` is \`false\`")
	@exit 1
else
	$(call green, "Found HOSTED_ZONE_ID: ${HOSTED_ZONE_ID}")
endif
	$(call blue, Validating Route53 configuration...)
	$(MAKE) env.validate.subdomain fqdn="${SUBDOMAIN}.${HOST_DOMAIN}."
	$(call green, Found configuration for ${SUBDOMAIN}.${HOST_DOMAIN})
ifeq ($(CREATE_VPC),false)
ifeq ($(VPC_ID),)
	$(call red, "\`VPC_ID\` must be set as an environment variable when \`CREATE_VPC\` is \`false\`")
	@exit 1
endif
ifeq ($(PUBLIC_SUBNET_ID),)
	$(call red, "\`PUBLIC_SUBNET_ID\` must be set as an environment variable when \`CREATE_VPC\` is \`false\`")
	@exit 1
endif
ifneq ($(PRIVATE_SUBNET_ID),)
	$(call red, "\`PRIVATE_SUBNET_ID\` must not be set when \`CREATE_VPC\` is \`false\`")
	@exit 1
endif
else ifeq ($(CREATE_VPC),true)
ifneq ($(VPC_ID),)
	$(call red, "\`VPC_ID\` must not be set as an environment variable when \`CREATE_VPC\` is \`true\`")
	@exit 1
endif
ifneq ($(PUBLIC_SUBNET_ID),)
	$(call red, "\`PUBLIC_SUBNET_ID\` must not be set as an environment variable when \`CREATE_VPC\` is \`true\`")
	@exit 1
endif
ifneq ($(PRIVATE_SUBNET_ID),)
	$(call red, "\`PRIVATE_SUBNET_ID\` must not be set when \`CREATE_VPC\` is \`false\`")
	@exit 1
endif
endif
endif

env.validate.boolean-vars:
	@$(foreach var,$(BOOLEAN_VARS),\
		if [ "$(value $(var))" != "" ] && [ "$(value $(var))" != "true" ] && [ "$(value $(var))" != "false" ]; then \
			echo "\033[0;31mERROR: \`$(var)\` must be unset, \`true\` or \`false\`\033[0m" && exit 1; \
		fi; \
	)

env.validate.vars:
	$(foreach var,$(REQUIRED_VARS),\
		$(if $(value $(var)),,$(error $(var) is not set. Please add $(var) to the environment variables.)))

env.validate.create-vpc.vars:
ifndef CREATE_VPC
	$(info 'CREATE_VPC' is not set. Defaulting to 'false')
	$(eval export CREATE_VPC := false)
	$(call blue, "**** This deployment uses an existing VPC ****")
endif
ifeq ($(CREATE_VPC),false)
	$(call blue, "**** This deployment uses an existing VPC ****")
ifeq ($(VPC_ID),)
	$(call red, "\`VPC_ID\` must be set as an environment variable when \`CREATE_VPC\` is \`false\`")
	@exit 1
else
	$(call green, "Found VPC_ID: ${VPC_ID}")
endif
ifeq ($(PUBLIC_SUBNET_ID),)
	$(call red, "\`PUBLIC_SUBNET_ID\` must be set as an environment variable when \`CREATE_VPC\` is \`false\`")
	@exit 1
else
	$(call green, "Found PUBLIC_SUBNET_ID: ${PUBLIC_SUBNET_ID}")
endif
else ifeq ($(CREATE_VPC),true)
	$(call blue, "**** This deployment includes a VPC ****")
endif

env.validate: check.dependencies env.validate.vars env.validate.boolean-vars env.validate.stage env.validate.create-vpc.vars env.validate.use-private-subnet.vars
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Found environment variables" 2>&1 | tee -a ${CFN_LOG_PATH}

options-screen:
	@echo "+-----------------------------------------------------------------------------------------+"
	@echo "|                                         \033[34mSuccess!\033[0m                                        |"
	@echo "+-----------------------------------------------------------------------------------------+"
	@echo "| \033[33mAvailable Actions:\033[0m                                                                      |"
	@echo "| * Run the pipeline: \`\033[96mSTAGE=<stage> make database.load.run releases=<releases>\033[0m\`          |"
	@echo "| * Load the database from backup: \`\033[96mSTAGE=<stage> make database.restore from_date=<date>\033[0m\` |"
	@echo "| * Log into the database: \`\033[96mSTAGE=<stage> make database.connect\033[0m\`                          |"
	@echo "| * Remove access services: \`\033[96mSTAGE=<stage> make infrastructure.access-services.delete\033[0m\`    |"
	@echo "+-----------------------------------------------------------------------------------------+"

infrastructure.deploy: 
	$(MAKE) -C ${APP_NAME}/infrastructure/ deploy

infrastructure.service.deploy:
	$(MAKE) -C ${APP_NAME}/infrastructure/ service.deploy

infrastructure.create-endpoint:
	$(MAKE) -C ${APP_NAME}/infrastructure/ service.deploy.create-endpoint service=$$service

infrastructure.access-services.deploy:
	$(MAKE) -C ${APP_NAME}/infrastructure/access-services/nat-gateway/ deploy
	$(MAKE) -C ${APP_NAME}/infrastructure/access-services/bastion-server/ deploy

infrastructure.access-services.nat-gateway.deploy:
	$(MAKE) -C ${APP_NAME}/infrastructure/access-services/nat-gateway/ deploy

infrastructure.access-services.bastion-server.deploy:
	$(MAKE) -C ${APP_NAME}/infrastructure/access-services/bastion-server/ deploy

infrastructure.access-services.bastion-server.connect:
	$(MAKE) -C ${APP_NAME}/infrastructure/access-services/bastion-server/ service.connect

database.deploy:
	$(MAKE) -C ${APP_NAME}/database/ deploy

database.service.deploy:
	$(MAKE) -C ${APP_NAME}/database/ service.deploy

database.connect:
	$(MAKE) infrastructure.access-services.bastion-server.connect

pipeline.deploy:
	$(MAKE) -C ${APP_NAME}/pipeline/ deploy

pipeline.service.deploy:
	$(MAKE) -C ${APP_NAME}/pipeline/ service.deploy

pipeline.jobs.deploy:
	$(MAKE) -C ${APP_NAME}/pipeline/ service.jobs.deploy

config.deploy:
	$(MAKE) -C ${APP_NAME}/pipeline/ service.config.deploy
	$(MAKE) -C ${APP_NAME}/database/ service.config.deploy

monitoring.create-subscriptions:
	$(MAKE) -C ${APP_NAME}/infrastructure service.monitoring.create-subscriptions

monitoring.subscribe-email:
	$(MAKE) -C ${APP_NAME}/infrastructure service.monitoring.subscribe-email

# TODO fix output & error handling
database.load.run: # args: align, kir, limit, releases
	@echo "Confirm payload:" && \
	[ "$$align" ] && align="$$align" || align=false && \
	[ "$$kir" ] && kir="$$kir" || kir=false && \
	[ "$$limit" ] && limit="$$limit" || limit="" && \
	[ "$$releases" ] && releases="$$releases" || releases="" && \
	[ "$$use_existing_build" ] && use_existing_build="$$use_existing_build" || use_existing_build=false && \
	[ "$$skip_load" ] && skip_load="$$skip_load" || skip_load=false && \
	payload="{\"align\":$$align,\"kir\":$$kir,\"limit\":\"$$limit\",\"releases\":\"$$releases\",\"mem_profile\":false,\"use_existing_build\":$$use_existing_build,\"skip_load\":$$skip_load}"&&\
	echo "$$payload" | jq -r && \
	echo "$$payload" | jq > payload.json
	@echo "Run pipeline with this payload? [y/N] \c " && read ans && [ $${ans:-N} = y ]
	@function_name="${STAGE}"-"${APP_NAME}"-"$$(cat ${FUNCTIONS_PATH}/environment.json | jq -r '.Functions.InvokePipeline.FunctionConfiguration.FunctionName')" && \
	echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Invoking $$function_name..." 2>&1 | tee -a ${CFN_LOG_PATH} && \
	echo "Payload:" >> ${CFN_LOG_PATH} && \
	cat payload.json >> ${CFN_LOG_PATH} && \
	aws lambda invoke \
		--cli-binary-format raw-in-base64-out \
		--function-name "$$function_name" \
		--payload file://payload.json \
		response.json \
		--output json  >> ${CFN_LOG_PATH} && \
	echo "Response:" && \
	echo "Response:" >> ${CFN_LOG_PATH} && \
	cat response.json | jq -r && \
	cat response.json | jq -r >> ${CFN_LOG_PATH} && \
	rm payload.json response.json

pipeline.invoke.validation-queries:
	@function_name="${STAGE}"-"${APP_NAME}"-"$$(cat ${FUNCTIONS_PATH}/environment.json | jq -r '.Functions.ExecuteValidationQueries.FunctionConfiguration.FunctionName')" && \
	echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Invoking $$function_name..." 2>&1 | tee -a ${CFN_LOG_PATH} && \
	aws lambda invoke \
		--cli-binary-format raw-in-base64-out \
		--function-name "$$function_name" \
		--log-type Tail \
		response.json \
		--output json  >> ${CFN_LOG_PATH} && \
	echo "Response:" >> ${CFN_LOG_PATH} && \
	cat response.json | jq -r && \
	cat response.json | jq -r >> ${CFN_LOG_PATH} && \
	rm response.json
	

# TODO database.load.status
# TODO database.load.abort
# TODO database.load.report

database.start:
	@echo "Starting $${APP_NAME} server..."
	@response=$$(aws ec2 start-instances --instance-ids ${INSTANCE_ID}) && \
	echo "Previous state: $$(echo "$$response" | jq -r '.StartingInstances[] | select(.InstanceId | contains("${INSTANCE_ID}")).PreviousState.Name')" && \
	echo "Current state: $$(echo "$$response" | jq -r '.StartingInstances[] | select(.InstanceId | contains("${INSTANCE_ID}")).CurrentState.Name')"

database.stop:
	@echo "Stopping $${APP_NAME} server..."
	@echo Instance ID: ${INSTANCE_ID}
	@response=$$(aws ec2 stop-instances --instance-ids ${INSTANCE_ID}) && \
	echo "Previous state: $$(echo "$$response" | jq -r '.StoppingInstances[] | select(.InstanceId | contains("${INSTANCE_ID}")).PreviousState.Name')" && \
	echo "Current state: $$(echo "$$response" | jq -r '.StoppingInstances[] | select(.InstanceId | contains("${INSTANCE_ID}")).CurrentState.Name')"

database.reboot:
	@echo "Rebooting $${APP_NAME} server..."
	@echo Instance ID: ${INSTANCE_ID}
	@response=$$(aws ec2 reboot-instances --instance-ids ${INSTANCE_ID}) && echo "$$response"
	$(MAKE) database.status

# TODO make sure database is running before syncing
database.sync-scripts:
	$(MAKE) -C ${APP_NAME}/database/ service.config.scripts.sync

# # TODO get expiration date, automate renewal
# database.ssl.get-expiration:
# 	$(MAKE) -C ${APP_NAME}/database/ service.ssl.get-expiration

database.config.update:
	@echo "Updating \`neo4j.conf\` up $${APP_NAME} server..."
	$(MAKE) -C ${APP_NAME}/database/ service.config.neo4j.update

database.ssl.renew-cert:
	$(MAKE) -C ${APP_NAME}/database/ service.ssl.renew-cert

database.backup:
	@echo "Backing up $${APP_NAME} server..."
	$(MAKE) -C ${APP_NAME}/database/ service.backup

database.backup.list:
	$(MAKE) -C ${APP_NAME}/database/ service.backup.list

# TODO call database.get.backups to list the available backups and prompt the user to select one
database.restore: #from_path=s3://<backup path>
	@echo "Restoring $${APP_NAME} data to server..."
	$(MAKE) -C ${APP_NAME}/database/ service.restore from_path=$$from_path

database.status:
	@aws ec2 describe-instances | \
		jq --arg iid "${INSTANCE_ID}" '.Reservations[].Instances[] | select(.InstanceId == $$iid) | {InstanceId, InstanceType, "Status": .State.Name, StateTransitionReason, ImageId}'

# TODO account for http or https and whether or not EIP or DNS is being used
database.get.endpoint:
	@echo "https://${SUBDOMAIN}.${HOST_DOMAIN}:7473/browser/"

database.get.credentials:
	@secret_string=$$(aws secretsmanager get-secret-value --secret-id ${APP_NAME}-${STAGE}-Neo4jCredentials | jq -r '.SecretString') && \
	echo "Username: $$(echo $$secret_string | jq -r '.NEO4J_USERNAME')" && \
	echo "Password: $$(echo $$secret_string | jq -r '.NEO4J_PASSWORD')"

database.get.private-ip:
	@private_ip=$$(aws ssm get-parameters \
		--names "/${APP_NAME}/${STAGE}/${AWS_REGION}/Neo4jPrivateIp" \
		--output json \
		| jq -r '.Parameters[0].Value') && \
	echo "$${private_ip}"

database.get.public-ip:
	@public_ip=$$(aws ssm get-parameters \
		--names "/${APP_NAME}/${STAGE}/${AWS_REGION}/Neo4jPublicIp" \
		--output json \
		| jq -r '.Parameters[0].Value') && \
	echo "$${public_ip}"

database.get.instance-id:
	@echo "${INSTANCE_ID}"

delete: # data=true/false ##=> Delete services
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting ${APP_NAME} in ${AWS_ACCOUNT}" 2>&1 | tee -a ${CFN_LOG_PATH}
	@[[ $$data != true ]] && echo "Data will not be deleted. To delete pass \`data=true\`" || true
	@echo "Delete all stacks from the \`${STAGE}\` environment? [y/N] \c " && read ans && [ $${ans:-N} = y ] && \
	if [ "${data}" = "true" ]; then \
		aws s3 rm --recursive s3://${DATA_BUCKET_NAME}; \
	fi
	$(MAKE) pipeline.delete
	$(MAKE) database.delete
	$(MAKE) infrastructure.delete
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Finished deleting ${APP_NAME} in ${AWS_ACCOUNT}" 2>&1 | tee -a ${CFN_LOG_PATH}


# Delete specific stacks
infrastructure.delete:
	$(MAKE) -C ${APP_NAME}/infrastructure/ delete

infrastructure.access-services.delete:
	$(MAKE) -C ${APP_NAME}/infrastructure/access-services/bastion-server/ delete
	$(MAKE) -C ${APP_NAME}/infrastructure/access-services/nat-gateway/ delete

infrastructure.access-services.bastion-server.delete:
	$(MAKE) -C ${APP_NAME}/infrastructure/access-services/bastion-server/ delete

infrastructure.access-services.nat-gateway.delete:
	$(MAKE) -C ${APP_NAME}/infrastructure/access-services/nat-gateway/ delete

infrastructure.delete-endpoint: #=> service=<string>
	$(MAKE) -C ${APP_NAME}/infrastructure/ service.delete.delete-endpoint service=$$service

database.delete:
	$(MAKE) -C ${APP_NAME}/database/ delete

pipeline.delete:
	$(MAKE) -C ${APP_NAME}/pipeline/ service.delete

pipeline.service.delete:
	$(MAKE) -C ${APP_NAME}/pipeline/ service.functions.delete

pipeline.jobs.delete:
	$(MAKE) -C ${APP_NAME}/pipeline/ service.jobs.delete

# Administrative functions
get.data: #=> Download the build data locally
	@mkdir -p ${ROOT_DIR}/data
	@aws s3 cp --recursive s3://${DATA_BUCKET_NAME}/data/ ${ROOT_DIR}/data/

get.logs: #=> Download all logs locally
	@aws s3 cp --recursive s3://${DATA_BUCKET_NAME}/logs/ ${LOGS_DIR}/

docs.build:
	@cd docs/ && make html

docs.url:
	@echo "$$(pwd)/docs/build/html/index.html"

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

	...::: Deploy the infrastructure :::...
	$ make infrastructure.deploy

	...::: Deploy the database :::...
	$ make database.deploy

	...::: Deploy the pipeline :::...
	$ make pipeline.deploy

	...::: Deploy the pipeline functions :::...
	$ make pipeline.service.deploy

	...::: Deploy the pipeline jobs :::...
	$ make pipeline.jobs.deploy

	...::: Deploy config files and scripts to S3 :::...
	$ make config.deploy

	...::: Run the StepFunctions State Machine to load Neo4j :::...
	$ make database.load releases=<version> align=<boolean> kir=<boolean> limit=<int>

	...::: Get Neo4j Credentials :::...
	$ make database.get.credentials

	...::: Get Neo4j Endpoint :::...
	$ make database.get.endpoint

	...::: Connect to the database server :::...
	$ make database.connect

	...::: Download CSV data from S3 :::...
	$ make get.data

	...::: Download logs from EC2 :::...
	$ make get.logs

	...::: Delete all CloudFormation based services and data :::...
	$ make delete

endef
