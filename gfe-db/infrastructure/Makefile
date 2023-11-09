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
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying ${SERVICE} service" 2>&1 | tee -a $${CFN_LOG_PATH}
	$(MAKE) service.deploy
ifeq ($(USE_PRIVATE_SUBNET),false)
	$(MAKE) service.deploy.update-dns
else ifeq ($(USE_PRIVATE_SUBNET),true)
ifeq ($(HAS_STAGE),null)
ifeq ($(CREATE_SSM_VPC_ENDPOINT),true)
	$(call blue, Creating Systems Manager VPC endpoint...)
	$(MAKE) service.deploy.create-endpoint service=ssm
endif
ifeq ($(CREATE_SECRETSMANAGER_VPC_ENDPOINT),true)
	$(call blue, Creating SecretsManager VPC endpoint...)
	$(MAKE) service.deploy.create-endpoint service=secretsmanager
endif
ifeq ($(CREATE_S3_VPC_ENDPOINT),true)
	$(call blue, Creating S3 VPC endpoint...)
	$(MAKE) service.deploy.create-endpoint service=s3
endif
endif
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

# TODO for CREATE_VPC=true and USE_PRIVATE_SUBNET=true, need to get the PrivateSubnetId from SSM Param Store, the endpoints are created
# but they are not associated with any subnet
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

# TODO test the return code of aws cloudformation deploy and exit if it is 1
service.deploy:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying CloudFormation" 2>&1 | tee -a $${CFN_LOG_PATH}
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
			AdminEmail="$${ADMIN_EMAIL}" \
			DataBucketName="$$DATA_BUCKET_NAME" \
			VpcId="$${VPC_ID}" \
			PublicSubnetId="$${PUBLIC_SUBNET_ID}" \
			PrivateSubnetId="$${PRIVATE_SUBNET_ID}" \
			GitHubPersonalAccessToken="$$GITHUB_PERSONAL_ACCESS_TOKEN"

service.nat-gateway.eip.create:
	@eip=$$(aws ssm get-parameters \
		--names "/${APP_NAME}/${STAGE}/${AWS_REGION}/NatGatewayElasticIp" \
		--output json | jq -r '.Parameters[0].Value') && \
	if [ "$$eip" = "null" ]; then \
		res=$$(aws ec2 allocate-address \
			--domain vpc \
			--tag-specifications "ResourceType=elastic-ip,Tags=[{Key=Name,Value=$${STAGE}-$${APP_NAME}-nat-gateway-eip}]" \
			) && \
		echo $$res | jq -r 2>&1 | tee -a $${CFN_LOG_PATH} && \
		eip=$$(echo $$res | jq -r '.PublicIp') && \
		echo "\033[0;34mSaving Elastic IP to SSM Parameter Store\033[0m" && \
		res=$$(aws ssm put-parameter \
			--name "/${APP_NAME}/${STAGE}/${AWS_REGION}/NatGatewayElasticIp" \
			--value "$$eip" \
			--type String \
			--overwrite) && \
		echo $$res | jq -r; \
	else \
		echo "\033[0;34mFound existing Elastic IP: $$eip\033[0m"; \
	fi

# wait for NAT Gateway to be created
service.nat-gateway.resource.create:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Creating NAT Gateway" 2>&1 | tee -a $${CFN_LOG_PATH}
	@nat_gateway_id=$$(aws ssm get-parameters \
		--names "/${APP_NAME}/${STAGE}/${AWS_REGION}/NatGatewayId" \
		--output json | jq -r '.Parameters[0].Value') && \
	if [ "$$nat_gateway_id" = "null" ]; then \
		eip=$$(aws ssm get-parameters \
			--names "/${APP_NAME}/${STAGE}/${AWS_REGION}/NatGatewayElasticIp" \
			--output json | jq -r '.Parameters[0].Value') && \
		if [ "$$eip" = "" ]; then \
			echo "\033[0;31mERROR: Elastic IP not found\033[0m" && exit 1; \
		else \
			allocation_id=$$(aws ec2 describe-addresses \
				--public-ips "$$eip" \
				--query "Addresses[0].AllocationId" \
				--output text) && \
			if [ "$$allocation_id" = "" ]; then \
				echo "\033[0;31mERROR: Allocation ID not found\033[0m" && exit 1; \
			fi; \
		fi && \
		public_subnet_id=$$(aws ssm get-parameters \
			--names "/${APP_NAME}/${STAGE}/${AWS_REGION}/PublicSubnetId" \
			--output json | jq -r '.Parameters[0].Value') && \
		if [ "$$public_subnet_id" = "null" ]; then \
			echo "\033[0;31mERROR: Public subnet id not found\033[0m" && exit 1; \
		fi && \
		echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying NAT Gateway" 2>&1 | tee -a $${CFN_LOG_PATH} && \
		res=$$(aws ec2 create-nat-gateway \
			--allocation-id "$$allocation_id" \
			--subnet-id "$$public_subnet_id" \
			--tag-specifications "ResourceType=natgateway,Tags=[{Key=Name,Value=$${STAGE}-$${APP_NAME}-nat-gateway}]") && \
		echo $$res | jq -r 2>&1 | tee -a $${CFN_LOG_PATH} && \
		nat_gateway_id=$$(echo $$res | jq -r '.NatGateway.NatGatewayId') && \
		res=$$(aws ec2 wait nat-gateway-available \
			--nat-gateway-ids "$$nat_gateway_id") && \
		echo "\033[0;34mSaving NAT Gateway ID to SSM Parameter Store\033[0m" && \
		res=$$(aws ssm put-parameter \
			--name "/${APP_NAME}/${STAGE}/${AWS_REGION}/NatGatewayId" \
			--value "$$nat_gateway_id" \
			--type String \
			--overwrite) && \
		echo $$res | jq -r; \
	else \
		echo "\033[0;34mFound existing NAT Gateway: $$nat_gateway_id\033[0m"; \
	fi

service.nat-gateway.deploy:
	$(MAKE) service.nat-gateway.eip.create
	$(MAKE) service.nat-gateway.resource.create


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
ifeq ($(CREATE_SSM_VPC_ENDPOINT),true)
	$(MAKE) service.delete.delete-endpoint service=ssm
endif
ifeq ($(CREATE_SECRETSMANAGER_VPC_ENDPOINT),true)
	$(MAKE) service.delete.delete-endpoint service=secretsmanager
endif
ifeq ($(CREATE_S3_VPC_ENDPOINT),true)
	$(MAKE) service.delete.delete-endpoint service=s3
endif
	$(MAKE) service.delete
	$(MAKE) service.parameters.delete

# TODO delete "/${APP_NAME}/${STAGE}/${AWS_REGION}/$${service}EndpointId"
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
	echo "Deleting VPC endpoint for $${service}..." && \
	res=$$(aws ec2 delete-vpc-endpoints \
		--vpc-endpoint-ids "$$endpoint_id") && \
	echo $$res | jq -r

service.delete:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting CloudFormation" 2>&1 | tee -a $${CFN_LOG_PATH}
	@aws cloudformation delete-stack \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" 2>&1 | tee -a $${CFN_LOG_PATH} || true && \
	aws cloudformation wait stack-delete-complete \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" 2>&1 | tee -a $${CFN_LOG_PATH} || true

service.nat-gateway.eip.delete:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting NAT Gateway Elastic IP" 2>&1 | tee -a $${CFN_LOG_PATH}
	@eip=$$(aws ssm get-parameters \
		--names "/${APP_NAME}/${STAGE}/${AWS_REGION}/NatGatewayElasticIp" \
		--output json | jq -r '.Parameters[0].Value') && \
	if [ "$$eip" != "null" ]; then \
		echo "\033[0;34mFound existing Elastic IP: $$eip\033[0m" && \
		allocation_id=$$(aws ec2 describe-addresses \
			--public-ips "$$eip" \
			--query "Addresses[0].AllocationId" \
			--output text) && \
		if [ "$$allocation_id" != "null" ]; then \
			res=$$(aws ec2 release-address \
				--allocation-id "$$allocation_id") && \
			echo $$res | jq -r 2>&1 | tee -a $${CFN_LOG_PATH} || true; \
		fi; \
		res=$$(aws ssm delete-parameter \
			--name "/${APP_NAME}/${STAGE}/${AWS_REGION}/NatGatewayElasticIp") && \
		echo $$res | jq -r && \
		echo "\033[0;34mSuccessfully deleted Elastic IP\033[0m"; \
	else \
		echo "\033[0;33mNo Elastic IP found\033[0m" && \
		exit 0; \
	fi

# wait for the NAT Gateway to be deleted before deleting the Elastic IP
service.nat-gateway.resource.delete:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting NAT Gateway" 2>&1 | tee -a $${CFN_LOG_PATH}
	@nat_gateway_id=$$(aws ssm get-parameters \
		--names "/${APP_NAME}/${STAGE}/${AWS_REGION}/NatGatewayId" \
		--output json | jq -r '.Parameters[0].Value') && \
		echo "$$nat_gateway_id" && \
	if [ "$$nat_gateway_id" != "null" ]; then \
		echo "\033[0;34mFound existing NAT Gateway: $$nat_gateway_id\033[0m" && \
		res=$$(aws ec2 delete-nat-gateway \
			--nat-gateway-id "$$nat_gateway_id") && \
		echo $$res | jq -r 2>&1 | tee -a $${CFN_LOG_PATH} || true; \
	else \
		echo "\033[0;33mNo NAT Gateway found\033[0m" && \
		exit 0; \
	fi && \
	echo "\033[0;34mWaiting for NAT Gateway to be deleted...\033[0m" && \
	res=$$(aws ec2 wait nat-gateway-deleted \
		--nat-gateway-ids "$$nat_gateway_id") && \
	echo $$res | jq -r 2>&1 | tee -a $${CFN_LOG_PATH} && \
	echo "\033[0;34mDeleting NAT Gateway ID from SSM Parameter Store\033[0m" && \
	res=$$(aws ssm delete-parameter \
		--name "/${APP_NAME}/${STAGE}/${AWS_REGION}/NatGatewayId") && \
	echo $$res | jq -r && \
	echo "\033[0;34mSuccessfully deleted NAT Gateway\033[0m"

service.nat-gateway.delete:
	$(MAKE) service.nat-gateway.resource.delete
	$(MAKE) service.nat-gateway.eip.delete

service.parameters.delete:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting SSM parameters" 2>&1 | tee -a ${CFN_LOG_PATH}
	@params=$$(aws ssm get-parameters-by-path \
		--path "/${APP_NAME}/${STAGE}/${AWS_REGION}" \
		--query "Parameters[*].Name" \
		--output json) && \
	if [ -n "$$params" ] && [ "$$params" != "[]" ]; then \
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