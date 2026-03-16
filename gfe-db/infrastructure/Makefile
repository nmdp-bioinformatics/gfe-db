.PHONY: service.update-vpc-endpoints
SERVICE := infrastructure

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

deploy: 
	@echo "\033[0;34m$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying ${SERVICE} services\033[0m" 2>&1 | tee -a $${CFN_LOG_PATH}
	$(MAKE) service.key-pair.create
	$(MAKE) service.deploy
	$(MAKE) service.access-services.deploy
	$(MAKE) service.monitoring.create-subscriptions

service.key-pair.create: ##=> Checks if the key pair already exists and creates it if it does not
	@echo "\033[0;34m$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Creating EC2 key pair\033[0m" 2>&1 | tee -a $${CFN_LOG_PATH}
	@key_pair="$$(aws ec2 describe-key-pairs --key-name ${EC2_KEY_PAIR_NAME} | jq '.KeyPairs[0].KeyName' || true)" && \
	[ "$$key_pair" ] && echo "Key pair found: $$key_pair" && exit 0 || echo "No key pair found..." && \
	echo "Creating EC2 key pair \"${EC2_KEY_PAIR_NAME}\"" && \
	aws ec2 create-key-pair --key-name ${EC2_KEY_PAIR_NAME} | jq -r '.KeyMaterial' > ${ROOT_DIR}/${EC2_KEY_PAIR_NAME}.pem && \
	echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Storing key pair ${EC2_KEY_PAIR_NAME} in Secrets Manager" 2>&1 | tee -a $$CFN_LOG_PATH && \
	aws secretsmanager create-secret --name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/EC2KeyPair" \
		--description "EC2 key pair for $${STAGE}-$${APP_NAME} database" \
		--secret-string "$$(cat ${ROOT_DIR}/${EC2_KEY_PAIR_NAME}.pem)" 2>&1 | tee -a $$CFN_LOG_PATH || true && \
	aws ssm put-parameter \
		--name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/EC2KeyPairName" \
		--type "String" \
		--value "${EC2_KEY_PAIR_NAME}" \
		--overwrite 2>&1 | tee -a $$CFN_LOG_PATH || true

service.deploy:
	@aws cloudformation deploy \
		--region $${AWS_REGION} \
		--template-file template.yaml \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" \
		--tags stage="$${STAGE}" app="$${APP_NAME}" service="${SERVICE}" branch="$$(git branch --show-current)" commit=$$(git rev-parse HEAD) \
		--parameter-overrides \
			Stage="$${STAGE}" \
			AppName="$${APP_NAME}" \
			createVpc="$${CREATE_VPC}" \
			deployNatGateway="$${DEPLOY_NAT_GATEWAY}" \
			AdminEmail="$${ADMIN_EMAIL}" \
			DataBucketName="$$DATA_BUCKET_NAME" \
			VpcId="$${VPC_ID}" \
			PublicSubnetId="$${PUBLIC_SUBNET_ID}" \
			PrivateSubnetId="$${PRIVATE_SUBNET_ID}" \
			ExternalNatGatewayId="$${EXTERNAL_NAT_GATEWAY_ID}" \
			ExternalPrivateRouteTableId="$${EXTERNAL_PRIVATE_ROUTE_TABLE_ID}" \
			GitHubPersonalAccessToken="$$GITHUB_PERSONAL_ACCESS_TOKEN"

service.access-services.deploy:
ifeq ($(DEPLOY_VPC_ENDPOINTS),true)
	$(MAKE) service.access-services.vpc-endpoints.deploy
else
	$(call blue, \`DEPLOY_VPC_ENDPOINTS\` is set to \`false\`. Update \`DEPLOY_VPC_ENDPOINTS\` to \`true\` to deploy.)
endif
ifeq ($(DEPLOY_NAT_GATEWAY),true)
	$(MAKE) service.access-services.nat-gateway.deploy
else
	$(call blue, \`DEPLOY_NAT_GATEWAY\` is set to \`false\`. Update \`DEPLOY_NAT_GATEWAY\` to \`true\` to deploy.)
endif
ifeq ($(DEPLOY_BASTION_SERVER),true)
	$(MAKE) service.access-services.bastion-server.deploy
else 
	$(call blue, \`DEPLOY_BASTION_SERVER\` is set to \`false\`. Update \`DEPLOY_BASTION_SERVER\` to \`true\` to deploy.)
endif

service.access-services.nat-gateway.deploy:
	$(MAKE) -C access-services/nat-gateway deploy

service.access-services.bastion-server.deploy:
	$(MAKE) -C access-services/bastion-server deploy

service.access-services.bastion-server.connect:
	$(MAKE) -C access-services/bastion-server service.bastion-server.connect

service.access-services.vpc-endpoints.deploy:
	$(MAKE) -C access-services/vpc-endpoints deploy

service.monitoring.create-subscriptions:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Creating SNS topic subscriptions" 2>&1 | tee -a $${CFN_LOG_PATH}
	@topic_arn=$$(aws ssm get-parameters \
		--names "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/DataPipelineErrorsTopicArn" \
		--with-decryption \
		--query "Parameters[0].Value" \
		--output text) && \
	for EMAIL in $$(echo $${SUBSCRIBE_EMAILS} | sed 's/,/ /g'); do \
		res=$$(aws sns subscribe \
			--topic-arn "$$topic_arn" \
			--protocol email \
			--notification-endpoint "$$EMAIL") && \
		echo $$res | jq -r 2>&1 | tee -a $${CFN_LOG_PATH} || true; \
	done

service.monitoring.subscribe-email: #email=<string>
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Creating SNS topic subscription" 2>&1 | tee -a $${CFN_LOG_PATH}
	@topic_arn=$$(aws ssm get-parameters \
		--names "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/DataPipelineErrorsTopicArn" \
		--with-decryption \
		--query "Parameters[0].Value" \
		--output text) && \
	aws sns subscribe \
		--topic-arn "$$topic_arn" \
		--protocol email \
		--notification-endpoint "$$email" 2>&1 | tee -a $${CFN_LOG_PATH} || true;

delete: ##=> Delete resources
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting ${SERVICE} service" 2>&1 | tee -a $${CFN_LOG_PATH}
ifeq ($(data),true)
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting data in $${DATA_BUCKET_NAME}" 2>&1 | tee -a $${CFN_LOG_PATH}
	@aws s3 rm --recursive s3://$${DATA_BUCKET_NAME}/
endif
ifeq ($(DEPLOY_BASTION_SERVER),true)
	$(MAKE) service.access-services.bastion-server.delete
endif
ifeq ($(DEPLOY_NAT_GATEWAY),true)
	$(MAKE) service.access-services.nat-gateway.delete
endif
ifeq ($(DEPLOY_VPC_ENDPOINTS),true)
	$(MAKE) service.access-services.vpc-endpoints.delete
endif
# endif
	$(MAKE) service.delete
	$(MAKE) service.key-pair.delete
	$(MAKE) service.parameters.delete

service.access-services.delete:
	$(MAKE) -C access-services/bastion-server delete
	$(MAKE) -C access-services/nat-gateway delete

service.access-services.bastion-server.delete:
	$(MAKE) -C access-services/bastion-server delete

service.access-services.nat-gateway.delete:
	$(MAKE) -C access-services/nat-gateway delete

service.access-services.vpc-endpoints.delete:
	$(MAKE) -C access-services/vpc-endpoints delete

service.delete:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting CloudFormation" 2>&1 | tee -a $${CFN_LOG_PATH}
	@aws cloudformation delete-stack \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" 2>&1 | tee -a $${CFN_LOG_PATH} || true && \
	aws cloudformation wait stack-delete-complete \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" 2>&1 | tee -a $${CFN_LOG_PATH} || true

service.key-pair.delete:
	@aws ec2 delete-key-pair --key-name "${EC2_KEY_PAIR_NAME}" 2>&1 | tee -a $$CFN_LOG_PATH || true && \
	aws ssm delete-parameter --name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/EC2KeyPairName" 2>&1 | tee -a $$CFN_LOG_PATH || true && \
	mv ${ROOT_DIR}/${EC2_KEY_PAIR_NAME}.pem ${ROOT_DIR}/deprecated-key-$${INSTANCE_ID}-$$(gdate -u +'%Y-%m-%d-%H-%M').pem || true && \
	echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting key pair ${EC2_KEY_PAIR_NAME} in Secrets Manager" 2>&1 | tee -a $$CFN_LOG_PATH && \
	aws secretsmanager delete-secret \
		--secret-id "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/EC2KeyPair" \
		--force-delete-without-recovery 2>&1 | tee -a $$CFN_LOG_PATH || true

service.parameters.delete:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting SSM parameters" 2>&1 | tee -a ${CFN_LOG_PATH}
	@params=$$(aws ssm get-parameters-by-path \
		--path "/${APP_NAME}/${STAGE}/${AWS_REGION}" \
		--query "Parameters[*].Name" \
		--output json | jq -r '.[]') && \
	echo $$params && \
	if [ -n "$$params" ] && [ "$$params" != "" ]; then \
		aws ssm delete-parameters --names $$params | jq -r; \
	else \
		echo "No parameters to delete."; \
	fi

#############
#  Helpers  #
#############

define HELP_MESSAGE

	Environment variables:

	SERVICE: "${SERVICE}"
		Description: Name of the service being deployed

	Common usage:

	...::: Deploy all CloudFormation based services :::...
	$ make deploy

	...::: Delete all CloudFormation based services :::...
	$ make delete

endef