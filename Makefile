# Environment variables
# include .env # Optional, include STAGE and AWS_PROFILE
include .env.${STAGE}
export

SPLASH_FONT := slant

export AWS_ACCOUNT = $(shell aws sts get-caller-identity \
	--query Account \
	--output text)
export HAS_STAGE = $(shell aws ssm get-parameters \
		--names "/${APP_NAME}/${STAGE}/${AWS_REGION}/Stage" \
		--output json \
		| jq -r '.Parameters[0].Value')

export ROOT_DIR = $(shell pwd)
export INFRA_DIR = ${ROOT_DIR}/${APP_NAME}/infrastructure
export DATABASE_DIR = ${ROOT_DIR}/${APP_NAME}/database
export PIPELINE_DIR = ${ROOT_DIR}/${APP_NAME}/pipeline
export LOGS_DIR = $(shell echo "${ROOT_DIR}/logs")
export CFN_LOG_PATH = $(shell echo "${LOGS_DIR}/cfn/logs.txt")
export PURGE_LOGS ?= false

# Conditionally required variable defaults
export CREATE_VPC ?= false
export USE_PRIVATE_SUBNET ?= false
export DEPLOY_NAT_GATEWAY ?=
export DEPLOY_BASTION_SERVER ?=
export CREATE_SSM_VPC_ENDPOINT ?=
export CREATE_SECRETSMANAGER_VPC_ENDPOINT ?=
export CREATE_S3_VPC_ENDPOINT ?=
export SSM_VPC_ENDPOINT_ID ?=
export SECRETSMANAGER_VPC_ENDPOINT_ID ?=
export S3_VPC_ENDPOINT_ID ?=
export VPC_ID ?=
export PUBLIC_SUBNET_ID ?=
export HOST_DOMAIN ?=
export SUBDOMAIN ?=
export PRIVATE_SUBNET_ID ?=
export ADMIN_IP ?=
export UBUNTU_AMI_ID ?= ami-0fc5d935ebf8bc3bc
export DATABASE_VOLUME_SIZE ?= 64

# Resource identifiers
export DATA_BUCKET_NAME ?= ${STAGE}-${APP_NAME}-${AWS_ACCOUNT}-${AWS_REGION}
export CONFIG_S3_PATH := config
export ECR_BASE_URI := ${AWS_ACCOUNT}.dkr.ecr.${AWS_REGION}.amazonaws.com
export BUILD_REPOSITORY_NAME ?= ${STAGE}-${APP_NAME}-build-service
export EC2_KEY_PAIR_NAME := $${STAGE}-$${APP_NAME}-$${AWS_REGION}-neo4j-key
export INSTANCE_ID = $(shell aws ssm get-parameters \
	--names "/${APP_NAME}/${STAGE}/${AWS_REGION}/Neo4jDatabaseInstanceId" \
	--output json \
	| jq -r '.Parameters[0].Value')

# Neo4j
export NEO4J_DATABASE_NAME ?= gfedb

# Required environment variables
REQUIRED_VARS := STAGE APP_NAME AWS_ACCOUNT AWS_REGION AWS_PROFILE SUBSCRIBE_EMAILS \
	GITHUB_REPOSITORY_OWNER GITHUB_REPOSITORY_NAME GITHUB_PERSONAL_ACCESS_TOKEN \
	ADMIN_EMAIL UBUNTU_AMI_ID NEO4J_PASSWORD GDS_VERSION

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

### Application Deployment Targets ###
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
	@awk -F'=' '{ \
		key = $$1; value = $$2; \
		if (key == "GITHUB_PERSONAL_ACCESS_TOKEN" || key == "DOCKER_PASSWORD" || key == "NEO4J_PASSWORD") value = "************"; \
		if (substr($$0, 1, 1) != "#") { \
			line = key "=" value; \
			line = substr(line, 1, 76); \
			if (length(line) > 76) line = line "..."; \
			printf "| %-79s |\n", line \
		} \
	}' .env.${STAGE}
	@echo "+---------------------------------------------------------------------------------+"
	@echo "\033[0;33mPlease confirm the above values are correct.\033[0m"

deploy: splash-screen logs.purge env.validate ##=> Deploy all services
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying ${APP_NAME} to ${AWS_ACCOUNT}" 2>&1 | tee -a ${CFN_LOG_PATH}
	$(MAKE) env.print
	@echo "Deploy stack to the \`${STAGE}\` environment? [y/N] \c " && read ans && [ $${ans:-N} = y ]
	$(MAKE) infrastructure.deploy 
	$(MAKE) database.deploy
	$(MAKE) pipeline.deploy
	@sh -c '$(MAKE) pipeline.state.build && $(MAKE) pipeline.state.load || echo "Pipeline state build failed"'
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Finished deploying ${APP_NAME}" 2>&1 | tee -a ${CFN_LOG_PATH}
	$(MAKE) options-screen

update: env.validate.stage env.validate
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Updating ${APP_NAME} to ${AWS_ACCOUNT}" 2>&1 | tee -a ${CFN_LOG_PATH}
	$(MAKE) env.print
	@echo "Update stack in the \`${STAGE}\` environment? [y/N] \c " && read ans && [ $${ans:-N} = y ]
	$(MAKE) infrastructure.deploy
	$(MAKE) database.deploy
	$(MAKE) pipeline.deploy
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Finished updating ${APP_NAME}" 2>&1 | tee -a ${CFN_LOG_PATH}

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

env.validate.use-private-subnet.vars:
ifeq ($(USE_PRIVATE_SUBNET),true)
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
ifeq ($(DEPLOY_NAT_GATEWAY),)
	$(call red, "\`DEPLOY_NAT_GATEWAY\` must be set when \`CREATE_VPC\` is \`false\` and \`USE_PRIVATE_SUBNET\` is \`true\`")
	@exit 1
else ifneq ($(DEPLOY_NAT_GATEWAY),)
	$(MAKE) env.validate.external-nat-gateway
endif
ifeq ($(CREATE_SSM_VPC_ENDPOINT),false)
ifeq ($(SSM_VPC_ENDPOINT_ID),)
	$(call red, "\`SSM_VPC_ENDPOINT_ID\` must be set as an environment variable when \`CREATE_VPC\` is \`true\` and \`CREATE_SSM_VPC_ENDPOINT\` is \`false\`")
	@exit 1
else
	$(call green, "Found SSM_VPC_ENDPOINT_ID: ${SSM_VPC_ENDPOINT_ID}")
endif
endif
ifeq ($(CREATE_SECRETSMANAGER_VPC_ENDPOINT),false)
ifeq ($(SECRETSMANAGER_VPC_ENDPOINT_ID),)
	$(call red, "\`SECRETSMANAGER_VPC_ENDPOINT_ID\` must be set as an environment variable when \`CREATE_VPC\` is \`true\` and \`CREATE_SSM_VPC_ENDPOINT\` is \`false\`")
	@exit 1
else
	$(call green, "Found SECRETSMANAGER_VPC_ENDPOINT_ID: ${SECRETSMANAGER_VPC_ENDPOINT_ID}")
endif
endif
ifeq ($(CREATE_S3_VPC_ENDPOINT),false)
ifeq ($(S3_VPC_ENDPOINT_ID),)
	$(call red, "\`S3_VPC_ENDPOINT_ID\` must be set as an environment variable when \`CREATE_VPC\` is \`true\` and \`CREATE_SSM_VPC_ENDPOINT\` is \`false\`")
	@exit 1
else
	$(call green, "Found S3_VPC_ENDPOINT_ID: ${S3_VPC_ENDPOINT_ID}")
endif
endif
else ifeq ($(CREATE_VPC),true)
ifneq ($(DEPLOY_NAT_GATEWAY),true)
	$(call red, "\`DEPLOY_NAT_GATEWAY\` must be set to \`true\` when \`CREATE_VPC\` is \`true\` and \`USE_PRIVATE_SUBNET\` is \`true\`")
	@exit 1
endif
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

env.validate.external-nat-gateway:
ifeq ($(DEPLOY_NAT_GATEWAY),false)
ifeq ($(EXTERNAL_NAT_GATEWAY_ID),)
	$(call red, "\`EXTERNAL_NAT_GATEWAY_ID\` must be set as an environment variable when \`DEPLOY_NAT_GATEWAY\` is \`false\`")
	@exit 1
else
	$(call green, "Found EXTERNAL_NAT_GATEWAY_ID: ${EXTERNAL_NAT_GATEWAY_ID}")
	$(call blue, Validating NAT Gateway configuration...)
	@res=$$(aws ec2 describe-route-tables \
		--filters "Name=vpc-id,Values=${VPC_ID}" | jq -r --arg SUBNET_ID "${PRIVATE_SUBNET_ID}" --arg NAT_ID "${EXTERNAL_NAT_GATEWAY_ID}" '.RouteTables[] \
			| select(.Associations[]? | .SubnetId == $$SUBNET_ID) | {RouteTableId: .RouteTableId, Routes: .Routes, IsNATGatewayRoutePresent: any(.Routes[]; .NatGatewayId == $$NAT_ID)}' \
	| jq -r '.IsNATGatewayRoutePresent') && \
	[[ $$res = "true" ]] && echo "\033[0;34mFound NAT Gateway route\033[0m" || (echo "\033[0;31mERROR: No NAT Gateway route found\033[0m" && exit 1)
endif
endif

env.validate.create-neo4j-users:
	@if [ -n "${CREATE_NEO4J_USERS}" ]; then \
	    valid_format=1; \
	    IFS=',' read -ra ADDR <<< "${CREATE_NEO4J_USERS}"; \
	    for user_pass in "$${ADDR[@]}"; do \
	        if [[ ! $$user_pass =~ ^[^:]+:[^:]+$$ ]]; then \
	            valid_format=0; \
	            break; \
	        fi; \
	    done; \
	    if [[ $$valid_format -eq 0 ]]; then \
	        echo "\033[0;31mERROR: Invalid Neo4j user format. Please use the format \`username:password\`.\033[0m"; \
	    fi; \
	else \
	    echo "\033[0;34mNo Neo4j users defined, skipping validation.\033[0m"; \
	fi

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

env.validate: check.dependencies env.validate.vars env.validate.boolean-vars env.validate.stage env.validate.create-vpc.vars env.validate.use-private-subnet.vars env.validate.create-neo4j-users
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Found environment variables" 2>&1 | tee -a ${CFN_LOG_PATH}

options-screen:
	@echo "+--------------------------------------------------------------------------------------------+"
	@echo "|                                         \033[34mSuccess!\033[0m                                           |"
	@echo "+--------------------------------------------------------------------------------------------+"
	@echo "| \033[33mAvailable Actions:\033[0m                                                                         |"
	@echo "| * Run the pipeline: \`\033[96mSTAGE=<stage> make database.load.run releases=<releases>\033[0m\`             |"
	@echo "| * Load the database from backup: \`\033[96mSTAGE=<stage> make database.restore from_path=<s3_path>\033[0m\` |"
	@echo "| * Log into the database: \`\033[96mSTAGE=<stage> make database.connect\033[0m\`                             |"
	@echo "| * Log into the Neo4j Browser: \`\033[96mSTAGE=<stage> make database.ui.connect\033[0m\`                     |"
	@echo "| * Remove access services: \`\033[96mSTAGE=<stage> make infrastructure.access-services.delete\033[0m\`       |"
	@echo "+--------------------------------------------------------------------------------------------+"

### Application Management Targets ###
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

monitoring.create-subscriptions:
	$(MAKE) -C ${APP_NAME}/infrastructure service.monitoring.create-subscriptions

monitoring.subscribe-email:
	$(MAKE) -C ${APP_NAME}/infrastructure service.monitoring.subscribe-email

database.deploy:
	$(MAKE) -C ${APP_NAME}/database/ deploy

database.service.deploy:
	$(MAKE) -C ${APP_NAME}/database/ service.deploy

database.config.deploy:
	$(MAKE) -C ${APP_NAME}/database/ service.config.deploy

database.connect:
ifeq ($(USE_PRIVATE_SUBNET),true)
	$(MAKE) infrastructure.access-services.bastion-server.connect
else
	$(MAKE) -C ${APP_NAME}/database/ service.connect
endif

database.ui.connect:
ifeq ($(USE_PRIVATE_SUBNET),true)
	$(MAKE) -C ${APP_NAME}/infrastructure/access-services/bastion-server/ service.ui.connect
else ifeq ($(USE_PRIVATE_SUBNET),false)
	$(MAKE) database.get.endpoint
endif

pipeline.deploy:
	$(MAKE) -C ${APP_NAME}/pipeline/ service.deploy

pipeline.service.update:
	$(MAKE) -C ${APP_NAME}/pipeline/ service.functions.deploy

pipeline.jobs.deploy:
	$(MAKE) -C ${APP_NAME}/pipeline/ service.jobs.deploy

pipeline.config.deploy:
	$(MAKE) -C ${APP_NAME}/pipeline/ service.config.deploy

pipeline.state.build:
	$(MAKE) -C ${APP_NAME}/pipeline/ service.state.build

pipeline.state.load: 
	$(MAKE) -C ${APP_NAME}/pipeline/ service.state.load

config.deploy:
	$(MAKE) database.config.deploy
	$(MAKE) pipeline.config.deploy

database.load.run: # args: align, kir, limit, releases
	@echo "Confirm payload:" && \
	[ "$$align" ] && align="$$align" || align=false && \
	[ "$$kir" ] && kir="$$kir" || kir=false && \
	[ "$$limit" ] && limit="$$limit" || limit=-1 && \
	[ "$$releases" ] && releases="$$releases" || releases="" && \
	[ "$$use_existing_build" ] && use_existing_build="$$use_existing_build" || use_existing_build=false && \
	[ "$$skip_load" ] && skip_load="$$skip_load" || skip_load=false && \
	payload="{\"align\":$$align,\"kir\":$$kir,\"limit\":\"$$limit\",\"releases\":\"$$releases\",\"mem_profile\":false,\"use_existing_build\":$$use_existing_build,\"skip_load\":$$skip_load}"&&\
	echo "$$payload" | jq -r && \
	echo "$$payload" | jq > payload.json
	@echo "Run pipeline with this payload? [y/N] \c " && read ans && [ $${ans:-N} = y ]
	@function_name="${STAGE}"-"${APP_NAME}"-"check-source-update" && \
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
	@function_name="${STAGE}"-"${APP_NAME}"-"$$(cat ${PIPELINE_DIR}/functions/environment.json | jq -r '.Functions.ExecuteValidationQueries.FunctionConfiguration.FunctionName')" && \
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

database.sync-scripts:
	$(MAKE) -C ${APP_NAME}/database/ service.config.scripts.sync

database.config.update:
	@echo "Updating \`neo4j.conf\` up $${APP_NAME} server..."
	$(MAKE) -C ${APP_NAME}/database/ service.config.neo4j.update

database.ssl.renew-cert:
	$(MAKE) -C ${APP_NAME}/database/ service.ssl.renew-cert

database.config.create-users:
	$(MAKE) -C ${APP_NAME}/database/ service.config.neo4j.create-users

database.backup:
	@echo "Backing up $${APP_NAME} server..."
	$(MAKE) -C ${APP_NAME}/database/ service.backup

database.get.current-backup:
	$(MAKE) -C ${APP_NAME}/database/ service.backup.get-current

database.backup.list:
	$(MAKE) -C ${APP_NAME}/database/ service.backup.list

database.restore: #from_path=s3://<backup path>
	@echo "Restoring $${APP_NAME} data to server..."
	$(MAKE) -C ${APP_NAME}/database/ service.restore from_path=$$from_path

database.status:
	@aws ec2 describe-instances | \
		jq --arg iid "${INSTANCE_ID}" '.Reservations[].Instances[] | (.InstanceId == $$iid) | {InstanceId, InstanceType, "Status": .State.Name, StateTransitionReason, ImageId}'

database.get.endpoint:
ifeq ($(USE_PRIVATE_SUBNET),true)
	@echo "http://localhost:7474/browser/"
else ifeq ($(USE_PRIVATE_SUBNET),false)
	@echo "https://${SUBDOMAIN}.${HOST_DOMAIN}:7473/browser/"
endif

database.get.credentials:
	@secret_string=$$(aws secretsmanager get-secret-value --secret-id /${APP_NAME}/${STAGE}/${AWS_REGION}/Neo4jCredentials | jq -r '.SecretString') && \
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

local.build:
	$(MAKE) -C ${APP_NAME}/local/ build

delete: # data=true/false ##=> Delete services
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting ${APP_NAME} in ${AWS_ACCOUNT}" 2>&1 | tee -a ${CFN_LOG_PATH}
	@[[ $$data != true ]] && echo "Data will not be deleted. To delete pass \`data=true\`" || true
	@echo "Delete all stacks from the \`${STAGE}\` environment? [y/N] \c " && read ans && [ $${ans:-N} = y ] && \
	if [ "${data}" = "true" ]; then \
		aws s3 rm --recursive --quiet s3://${DATA_BUCKET_NAME}; \
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

	************************
	* Deployment Variables *
	************************

	STAGE: "${STAGE}"
		Description: (string) Feature branch name used as part of stacks name

	AWS_PROFILE: "${AWS_PROFILE}"
		Description: (string) AWS profile to use for deployment

	APP_NAME: "${APP_NAME}"
		Description: (string) Stack Name already deployed

	AWS_REGION: "${AWS_REGION}":
		Description: (string) AWS region for deployment

	CREATE_VPC: "${CREATE_VPC}"
		Description: (boolean) Create a new VPC or use an existing one

	USE_PRIVATE_SUBNET: "${USE_PRIVATE_SUBNET}"
		Description: (boolean) Use a private subnet for Neo4j

	CREATE_SSM_VPC_ENDPOINT: "${CREATE_SSM_VPC_ENDPOINT}"
		Description: (boolean) Create a new SSM VPC Endpoint or use an existing one

	SSM_VPC_ENDPOINT_ID: "${SSM_VPC_ENDPOINT_ID}"
		Description: (string) ID of an existing SSM VPC Endpoint, required when CREATE_SSM_VPC_ENDPOINT is false
		and USE_PRIVATE_SUBNET is true

	CREATE_SECRETSMANAGER_VPC_ENDPOINT: "${CREATE_SECRETSMANAGER_VPC_ENDPOINT}"
		Description: (boolean) Create a new SecretsManager VPC Endpoint or use an existing one

	SECRETSMANAGER_VPC_ENDPOINT_ID: "${SECRETSMANAGER_VPC_ENDPOINT_ID}"
		Description: (string) ID of an existing SecretsManager VPC Endpoint, required when CREATE_SECRETSMANAGER_VPC_ENDPOINT is false
		and USE_PRIVATE_SUBNET is true

	CREATE_S3_VPC_ENDPOINT: "${CREATE_S3_VPC_ENDPOINT}"
		Description: (boolean) Create a new S3 VPC Endpoint or use an existing one, required when USE_PRIVATE_SUBNET is true

	S3_VPC_ENDPOINT_ID: "${S3_VPC_ENDPOINT_ID}"
		Description: (string) ID of an existing S3 VPC Endpoint, required when CREATE_S3_VPC_ENDPOINT is false
		and USE_PRIVATE_SUBNET is true

	DEPLOY_NAT_GATEWAY: "${DEPLOY_NAT_GATEWAY}"
		Description: (boolean) Deploy a NAT Gateway or use an existing one, required when USE_PRIVATE_SUBNET is true

	EXTERNAL_NAT_GATEWAY_ID: "${EXTERNAL_NAT_GATEWAY_ID}"
		Description: (string) ID of an existing NAT Gateway, required when DEPLOY_NAT_GATEWAY is false
		and USE_PRIVATE_SUBNET is true

	DEPLOY_BASTION_SERVER: "${DEPLOY_BASTION_SERVER}"
		Description: (boolean) Deploy a Bastion Server or use an existing one, required when USE_PRIVATE_SUBNET is true
		
	ADMIN_IP: "${ADMIN_IP}"
		Description: (string) IP address to allow SSH access to the Bastion Server, required when DEPLOY_BASTION_SERVER is true

	ADMIN_EMAIL: "${ADMIN_EMAIL}"
		Description: (string) Admin email address for Neo4j server SSL certificate management, required when USE_PRIVATE_SUBNET is false

	SUBSCRIBE_EMAILS: "${SUBSCRIBE_EMAILS}"
		Description: (string) Comma separated list of email addresses to subscribe to CloudWatch notifications

	HOST_DOMAIN: "${HOST_DOMAIN}"
		Description: (string) Domain name for the Neo4j server, required when USE_PRIVATE_SUBNET is false

	SUBDOMAIN: "${SUBDOMAIN}"
		Description: (string) Subdomain name for the Neo4j server, required when USE_PRIVATE_SUBNET is false

	HOSTED_ZONE_ID: "${HOSTED_ZONE_ID}"
		Description: (string) Route53 hosted zone ID, required when USE_PRIVATE_SUBNET is false

	VPC_ID: "${VPC_ID}"
		Description: (string) ID of an existing VPC, required when CREATE_VPC is false

	PUBLIC_SUBNET_ID: "${PUBLIC_SUBNET_ID}"
		Description: (string) ID of an existing public subnet, required when CREATE_VPC is false

	PRIVATE_SUBNET_ID: "${PRIVATE_SUBNET_ID}"
		Description: (string) ID of an existing private subnet, required when CREATE_VPC is false
		and USE_PRIVATE_SUBNET is true

	UBUNTU_AMI_ID: "${UBUNTU_AMI_ID}"
		Description: (string) ID of an existing AMI for Ubuntu 22.04

	APOC_VERSION: "${APOC_VERSION}"
		Description: (string) Version of APOC to install

	GDS_VERSION: "${GDS_VERSION}"
		Description: (string) Version of Neo4j Graph Data Science plugin to install

	GITHUB_PERSONAL_ACCESS_TOKEN: "${GITHUB_PERSONAL_ACCESS_TOKEN}"
		Description: (string) GitHub personal access token for downloading releases

	GITHUB_REPOSITORY_OWNER: "${GITHUB_REPOSITORY_OWNER}"
		Description: (string) GitHub repository owner for downloading releases

	GITHUB_REPOSITORY_NAME: "${GITHUB_REPOSITORY_NAME}"
		Description: (string) GitHub repository name for downloading releases

	FEATURE_SERVICE_URL: "${FEATURE_SERVICE_URL}"
		Description: (string) URL of the Feature Service API

	AWS_ACCOUNT: "${AWS_ACCOUNT}":
		Description: AWS account ID for deployment

	*************************
	* Application variables *
	*************************

	AWS_ACCOUNT: "${AWS_ACCOUNT}"
		Description: AWS account ID for deployment

	HAS_STAGE: "${HAS_STAGE}"
		Description: (boolean) Whether or not the stack has a stage already deployed

	ROOT_DIR: "${ROOT_DIR}"
		Description: Root directory of the project

	INFRA_DIR: "${INFRA_DIR}"
		Description: Path to the infrastructure directory

	DATABASE_DIR: "${DATABASE_DIR}"
		Description: Path to the database directory

	PIPELINE_DIR: "${PIPELINE_DIR}"
		Description: Path to the pipeline directory

	DATA_BUCKET_NAME "${DATA_BUCKET_NAME}"
		Description: Name of the S3 bucket for data, config, logs and backups

	ECR_BASE_URI: "${ECR_BASE_URI}"
		Description: Base URI of the ECR repository for the build job Application

	BUILD_REPOSITORY_NAME: "${BUILD_REPOSITORY_NAME}"
		Description: Name of the ECR repository for the build job Application

	EC2_KEY_PAIR_NAME: "${EC2_KEY_PAIR_NAME}"
		Description: Name of the EC2 key pair for the Bastion and Neo4j servers

	INSTANCE_ID: "${INSTANCE_ID}"
		Description: ID of the Neo4j EC2 instance

	PIPELINE_STATE_PATH: "${PIPELINE_STATE_PATH}"
		Description: S3 path to the pipeline state file

	FUNCTIONS_PATH: "${FUNCTIONS_PATH}"
		Description: Path to the Lambda functions directory

	****************
	* Common Usage *
	****************

	...::: Deploy all CloudFormation based services :::...
	$ make deploy

	...::: Delete all CloudFormation based services including data :::...
	$ make delete data=true

	...::: Deploy the infrastructure layer :::...
	$ make infrastructure.deploy

	...::: Delete the infrastructure layer :::...
	$ make infrastructure.delete

	...::: Update only the infrastructure CloudFormation :::...
	$ make infrastructure.service.deploy

	...::: Create an endpoint :::...
	$ make infrastructure.create-endpoint service=<ssm|secretmanager|s3>

	...::: Delete an endpoint :::...
	$ make infrastructure.delete-endpoint service=<ssm|secretmanager|s3>

	...::: Deploy all access services :::...
	$ make infrastructure.access-services.deploy

	...::: Delete all access services :::...
	$ make infrastructure.access-services.delete

	...::: Deploy the NAT Gateway access service :::...
	$ make infrastructure.access-services.nat-gateway.deploy

	...::: Delete the NAT Gateway access service :::...
	$ make infrastructure.access-services.nat-gateway.delete

	...::: Deploy the Bastion Server access service :::...
	$ make infrastructure.access-services.bastion-server.deploy

	...::: Delete the Bastion Server access service :::...
	$ make infrastructure.access-services.bastion-server.delete

	...::: Create CloudWatch subscriptions :::...
	$ make monitoring.create-subscriptions

	...::: Subscribe an email to CloudWatch notifications :::...
	$ make monitoring.subscribe-email

	...::: Connect to the Bastion Server :::...
	$ make infrastructure.access-services.bastion-server.connect

	...::: Deploy the database layer :::...
	$ make database.deploy

	...::: Delete the database layer :::...
	$ make database.delete

	...::: Update only the database CloudFormation :::...
	$ make database.service.deploy

	...::: Connect to the database instance :::...
	$ make database.connect

	...::: Connect to the Neo4j Browser :::...
	$ make database.ui.connect

	...::: Run the StepFunctions State Machine to load Neo4j :::...
	$ make database.load releases=<version> align=<boolean> kir=<boolean> limit=<int>

	...::: Get Neo4j Credentials :::...
	$ make database.get.credentials

	...::: Get Neo4j Endpoint :::...
	$ make database.get.endpoint

	...::: Get Neo4j Private IP :::...
	$ make database.get.private-ip

	...::: Get Neo4j Instance ID :::...
	$ make database.get.instance-id

	...::: Stop the database instance :::...
	$ make database.stop

	...::: Start the database instance :::...
	$ make database.start

	...::: Reboot the database instance :::...
	$ make database.reboot

	...::: Sync bash scripts to the database instance :::...
	$ make database.sync-scripts

	...::: Update the Neo4j configuration :::...
	$ make database.config.update

	...::: Renew the Neo4j SSL certificate :::...
	$ make database.ssl.renew-cert

	...::: Backup the database :::...
	$ make database.backup

	...::: List database backups :::...
	$ make database.backup.list

	...::: Restore the database from a backup :::...
	$ make database.restore from_path=<s3_path>

	...::: Describe the database status :::...
	$ make database.status

	...::: Deploy the pipeline :::...
	$ make pipeline.deploy

	...::: Delete the pipeline :::...
	$ make pipeline.delete

	...::: Update only the pipeline CloudFormation including Lambda functions :::...
	$ make pipeline.service.update

	...::: Deploy the pipeline jobs as Docker images to ECR:::...
	$ make pipeline.jobs.deploy

	...::: Delete the pipeline jobs :::...
	$ make pipeline.jobs.delete

	...::: Deploy all config files and scripts to S3 :::...
	$ make config.deploy

	...::: Deploy only the database config files and scripts to S3 :::...
	$ make database.config.deploy

	...::: Download CSV data from S3 :::...
	$ make get.data

	...::: Download server logs from the Neo4j instance :::...
	$ make get.logs

	...::: Build the documentation :::...
	$ make docs.build

	...::: Get the documentation URL :::...
	$ make docs.url

endef
