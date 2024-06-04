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
ifeq ($(USE_PRIVATE_SUBNET),true)
	$(MAKE) service.access-services.deploy
else ifeq ($(USE_PRIVATE_SUBNET),false)
	$(MAKE) service.deploy.update-dns
# else ifeq ($(USE_PRIVATE_SUBNET),true)
# ifeq ($(HAS_STAGE),null)
# ifeq ($(SKIP_CONFIGURE_VPC_ENDPOINTS),false)
# ifeq ($(CREATE_SSM_VPC_ENDPOINT),true)
# 	$(call blue, Creating Systems Manager VPC endpoint...)
# 	$(MAKE) service.deploy.create-endpoint service=ssm
# else 
# 	$(MAKE) service.deploy.create-endpoint-id-parameter service=ssm vpce_id=${SSM_VPC_ENDPOINT_ID}
# 	@vpce_type=$$(aws ec2 describe-vpc-endpoints --query "VpcEndpoints[?VpcEndpointId=='${SSM_VPC_ENDPOINT_ID}'].VpcEndpointType" | jq -r '.[]') && \
# 	echo "vpce_type: $$vpce_type" && \
# 	if [ "$$vpce_type" = "Interface" ]; then \
# 		echo "\033[0;34mAdding Neo4j Security Group to existing Systems Manager VPC Interface endpoint...\033[0m" && \
# 		$(MAKE) service.deploy.modify-endpoint service=ssm; \
# 	fi
# endif
# ifeq ($(CREATE_SECRETSMANAGER_VPC_ENDPOINT),true)
# 	$(call blue, Creating SecretsManager VPC endpoint...)
# 	$(MAKE) service.deploy.create-endpoint service=secretsmanager
# else 
# 	$(MAKE) service.deploy.create-endpoint-id-parameter service=secretsmanager vpce_id=${SECRETSMANAGER_VPC_ENDPOINT_ID}
# 	@vpce_type=$$(aws ec2 describe-vpc-endpoints --query "VpcEndpoints[?VpcEndpointId=='${SECRETSMANAGER_VPC_ENDPOINT_ID}'].VpcEndpointType" | jq -r '.[]') && \
# 	echo "vpce_type: $$vpce_type" && \
# 	if [ "$$vpce_type" = "Interface" ]; then \
# 		echo "\033[0;34mAdding Neo4j Security Group to existing SecretsManager VPC Interface endpoint...\033[0m" && \
# 		$(MAKE) service.deploy.modify-endpoint service=secretsmanager; \
# 	fi
# endif
# ifeq ($(CREATE_S3_VPC_ENDPOINT),true)
# 	$(call blue, Creating S3 VPC endpoint...)
# 	$(MAKE) service.deploy.create-endpoint service=s3
# else
# 	$(MAKE) service.deploy.create-endpoint-id-parameter service=s3 vpce_id=${S3_VPC_ENDPOINT_ID}
# 	@vpce_type=$$(aws ec2 describe-vpc-endpoints --query "VpcEndpoints[?VpcEndpointId=='${S3_VPC_ENDPOINT_ID}'].VpcEndpointType" | jq -r '.[]') && \
# 	echo "vpce_type: $$vpce_type" && \
# 	if [ "$$vpce_type" = "Interface" ]; then \
# 		echo "\033[0;34mAdding Neo4j Security Group to existing S3 VPC Interface endpoint...\033[0m" && \
# 		$(MAKE) service.deploy.modify-endpoint service=s3; \
# 	fi
# endif
# else
# 	$(call blue, Skipping VPC endpoint configuration...)
# endif
# endif
endif
	$(MAKE) service.monitoring.create-subscriptions

service.deploy.update-dns:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Updating DNS records" 2>&1 | tee -a $${CFN_LOG_PATH}
	@config_path=${INFRA_DIR}/change-batch.json && \
	elastic_ip=$$(aws ssm get-parameters \
		--names "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/Neo4jDatabaseEndpoint" \
		--with-decryption \
		--query "Parameters[0].Value" \
		--output text) && \
	sed -e "s/<domain_name>/${SUBDOMAIN}.${HOST_DOMAIN}./g" -e "s/<elastic_ip>/$$elastic_ip/g" $$config_path > $$config_path.tmp && \
	echo "Updating DNS records with:" && \
	cat $$config_path.tmp | jq -r && \
	res=$$(aws route53 change-resource-record-sets --hosted-zone-id $${HOSTED_ZONE_ID} --change-batch file://$$config_path.tmp) && \
	echo && \
	echo "Waiting for update to complete:" && \
	echo $$res | jq -r && \
	aws route53 wait resource-record-sets-changed --id $$(echo $$res | jq -r '.ChangeInfo.Id') && \
	echo && \
	echo "Successfully updated DNS records"

service.deploy.create-endpoint: #=> service=<string>
	@if [ -z $$service ]; then \
		echo "\033[0;31mERROR: No service specified\033[0m" && exit 1; \
	elif [ $$service != "ssm" ] && [ $$service != "secretsmanager" ] && [ $$service != "s3" ]; then \
		echo "\033[0;31mERROR: Invalid service specified. Valid services are 'ssm', 'secretsmanager' or 's3'\033[0m" && exit 1; \
	fi && \
	vpc_id=$$(aws ssm get-parameters \
		--names "/${APP_NAME}/${STAGE}/${AWS_REGION}/VpcId" \
		--output json | jq -r '.Parameters[0].Value') && \
	if [ "$$vpc_id" = "" ]; then \
		echo "\033[0;31mERROR: VPC ID not found\033[0m" && exit 1; \
	fi && \
	neo4j_sg_id=$$(aws ssm get-parameters \
		--names "/${APP_NAME}/${STAGE}/${AWS_REGION}/Neo4jDatabaseSecurityGroupId" \
        --output json | jq -r '.Parameters[0].Value') && \
	if [ "$$neo4j_sg_id" = "" ]; then \
		echo "\033[0;31mERROR: Neo4j security group id not found\033[0m" && exit 1; \
	fi && \
	private_subnet_id=$$(aws ssm get-parameters \
		--names "/${APP_NAME}/${STAGE}/${AWS_REGION}/PrivateSubnetId" \
		--output json | jq -r '.Parameters[0].Value') && \
	if [ "$$private_subnet_id" = "" ]; then \
		echo "\033[0;31mERROR: Private subnet id not found\033[0m" && exit 1; \
	fi && \
	if [ $$service = "ssm" ] && [ ${CREATE_SSM_VPC_ENDPOINT} != "true" ]; then \
		echo "\033[0;31mERROR: \`CREATE_SSM_VPC_ENDPOINT\` must be set to \`true\` to create a Systems Manager VPC endpoint\033[0m" && exit 1; \
	fi && \
	if [ $$service = "secretsmanager" ] && [ ${CREATE_SECRETSMANAGER_VPC_ENDPOINT} != "true" ]; then \
		echo "\033[0;31mERROR: \`CREATE_SECRETSMANAGER_VPC_ENDPOINT\` must be set to \`true\` to create a SecretsManager VPC endpoint\033[0m" && exit 1; \
	fi && \
	if [ $$service = "s3" ] && [ ${CREATE_S3_VPC_ENDPOINT} != "true" ]; then \
		echo "\033[0;31mERROR: \`CREATE_S3_VPC_ENDPOINT\` must be set to \`true\` to create an S3 VPC endpoint\033[0m" && exit 1; \
	fi && \
	res=$$(aws ec2 create-vpc-endpoint \
		--vpc-id $$vpc_id \
		--service-name com.amazonaws.$${AWS_REGION}.$$service \
		--subnet-ids $$private_subnet_id \
		--security-group-ids $$neo4j_sg_id \
		--vpc-endpoint-type Interface \
		--tag-specifications "ResourceType=vpc-endpoint,Tags=[{Key=Name,Value=$${STAGE}-$${APP_NAME}-$${service}-vpce}]") && \
	vpc_endpoint_id=$$(echo $$res | jq -r '.VpcEndpoint.VpcEndpointId') && \
	echo "Saving endpoint id to SSM Parameter Store" && \
	res=$$(aws ssm put-parameter \
		--name "/${APP_NAME}/${STAGE}/${AWS_REGION}/$${service}EndpointId" \
		--value "$$vpc_endpoint_id" \
		--type String \
		--overwrite) && \
	echo $$res | jq -r
	$(call blue, Successfully created VPC endpoint for $$service)

service.deploy.create-endpoint-id-parameter: #=> service=<service> vpce_id=<string>
	$(call blue, Updating SSM Parameters for $$service VPC endpoint...)
	@res=$$(aws ssm put-parameter \
		--name /${APP_NAME}/${STAGE}/${AWS_REGION}/"$$service"EndpointId \
		--value "$$vpce_id" \
		--type String \
		--overwrite) && \
	echo $$res | jq -r

service.deploy.modify-endpoint: #=> service=<service>
	@if [ -z $$service ]; then \
		echo "\033[0;31mERROR: No service specified\033[0m" && exit 1; \
	elif [ $$service != "ssm" ] && [ $$service != "secretsmanager" ] && [ $$service != "s3" ]; then \
		echo "\033[0;31mERROR: Invalid service specified. Valid services are 'ssm', 'secretsmanager' or 's3'\033[0m" && exit 1; \
	fi && \
	neo4j_sg_id=$$(aws ssm get-parameters \
		--names "/${APP_NAME}/${STAGE}/${AWS_REGION}/Neo4jDatabaseSecurityGroupId" \
        --output json | jq -r '.Parameters[0].Value') && \
	if [ "$$neo4j_sg_id" = "" ]; then \
		echo "\033[0;31mERROR: Neo4j security group id not found\033[0m" && exit 1; \
	fi && \
	res=$$(aws ec2 modify-vpc-endpoint \
		--vpc-endpoint-id $$(aws ssm get-parameters \
			--names "/${APP_NAME}/${STAGE}/${AWS_REGION}/$${service}EndpointId" \
			--output json | jq -r '.Parameters[0].Value') \
		--add-security-group-ids $$neo4j_sg_id) && \
	echo $$res | jq -r
	$(call blue, Successfully modified VPC endpoint for $$service)

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
			usePrivateSubnet="$${USE_PRIVATE_SUBNET}" \
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
ifeq ($(USE_PRIVATE_SUBNET),true)
ifeq ($(DEPLOY_NAT_GATEWAY),true)
	$(MAKE) service.access-services.nat-gateway.delete
endif
ifeq ($(DEPLOY_VPC_ENDPOINTS),true)
	$(MAKE) service.access-services.vpc-endpoints.delete
endif
ifeq ($(DEPLOY_BASTION_SERVER),true)
	$(MAKE) service.access-services.bastion-server.delete
endif
# ifeq ($(CREATE_SSM_VPC_ENDPOINT),true)
# 	$(MAKE) service.delete.delete-endpoint service=ssm
# else
# 	$(MAKE) service.delete.modify-endpoint service=ssm
# endif
# ifeq ($(CREATE_SECRETSMANAGER_VPC_ENDPOINT),true)
# 	$(MAKE) service.delete.delete-endpoint service=secretsmanager
# else
# 	$(MAKE) service.delete.modify-endpoint service=secretsmanager
# endif
# ifeq ($(CREATE_S3_VPC_ENDPOINT),true)
# 	$(MAKE) service.delete.delete-endpoint service=s3
# else
# 	$(MAKE) service.delete.modify-endpoint service=s3
# endif
endif
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

service.delete.delete-endpoint: #=> service=<string>
	@if [ -z $$service ]; then \
		echo "No service specified" && exit 1; \
	elif [ $$service != "ssm" ] && [ $$service != "secretsmanager" ] && [ $$service != "s3" ]; then \
		echo "\033[0;31mERROR: Invalid service specified. Valid services are 'ssm', 'secretsmanager' or 's3'\033[0m" && exit 1; \
	fi
	@-endpoint_id=$$(aws ssm get-parameters \
		--names "/${APP_NAME}/${STAGE}/${AWS_REGION}/$${service}EndpointId" \
		--with-decryption \
		--query "Parameters[0].Value" \
		--output text) && \
	if [ "$$endpoint_id" != "null" ] && [ "$$endpoint_id" != "None" ]; then \
		echo "Deleting VPC endpoint for $${service}..." && \
		res=$$(aws ec2 delete-vpc-endpoints \
			--vpc-endpoint-ids "$$endpoint_id") && \
		echo $$res | jq -r && \
		echo "Deleting SSM parameter..." && \
		res=$$(aws ssm delete-parameter \
			--name "/${APP_NAME}/${STAGE}/${AWS_REGION}/$${service}EndpointId") && \
		echo $$res | jq -r; \
	else \
		echo "\033[0;33mNo endpoint ID found for $${service}\033[0m"; \
	fi

service.delete.modify-endpoint: #=> service=<service>
	@-if [ -z $$service ]; then \
		echo "\033[0;31mERROR: No service specified\033[0m" && exit 1; \
	elif [ $$service != "ssm" ] && [ $$service != "secretsmanager" ] && [ $$service != "s3" ]; then \
		echo "\033[0;31mERROR: Invalid service specified. Valid services are 'ssm', 'secretsmanager' or 's3'\033[0m" && exit 1; \
	fi && \
	res=$$(aws ec2 modify-vpc-endpoint \
		--vpc-endpoint-id $$(aws ssm get-parameters \
			--names "/${APP_NAME}/${STAGE}/${AWS_REGION}/$${service}EndpointId" \
			--output json | jq -r '.Parameters[0].Value') \
		--remove-security-group-ids $$(aws ssm get-parameters \
			--names "/${APP_NAME}/${STAGE}/${AWS_REGION}/Neo4jDatabaseSecurityGroupId" \
			--output json | jq -r '.Parameters[0].Value')) && \
	echo $$res | jq -r && \
	echo "Successfully removed Neo4j security group from $${service} VPC endpoint"

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