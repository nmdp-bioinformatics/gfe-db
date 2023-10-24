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
endif
	$(MAKE) service.monitoring.create-subscriptions
ifeq ($(HAS_STAGE),null)
ifeq ($(CREATE_SSM_VPC_ENDPOINT),true)
	$(call blue, Creating Systems Manager VPC endpoint..)
	$(MAKE) service.deploy.create-endpoint service=ssm
endif
ifeq ($(CREATE_SECRETSMANAGER_VPC_ENDPOINT),true)
	$(call blue, Creating SecretsManager VPC endpoint...)
	$(MAKE) service.deploy.create-endpoint service=secretsmanager
endif
endif

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
	@neo4j_sg_id=$$(aws ssm get-parameters \
		--names "/${APP_NAME}/${STAGE}/${AWS_REGION}/Neo4jDatabaseSecurityGroupId" \
        --output json | jq -r '.Parameters[0].Value') && \
	if [ "$$neo4j_sg_id" = "" ]; then \
		echo "Neo4j security group id not found" && exit 1; \
	fi && \
	if [ -z $$service ]; then \
		echo "No service specified" && exit 1; \
	elif [ $$service != "ssm" ] && [ $$service != "secretsmanager" ]; then \
		echo "Invalid service specified. Valid services are 'ssm', 'secretsmanager'" && exit 1; \
	fi && \
	res=$$(aws ec2 create-vpc-endpoint \
		--vpc-id $${VPC_ID} \
		--service-name com.amazonaws.$${AWS_REGION}.$$service \
		--security-group-ids $$neo4j_sg_id \
		--vpc-endpoint-type Interface) \
		--tag-specification "ResourceType=vpc-endpoint,Tags=[{Key=Name,Value=$${STAGE}-$${APP_NAME}-$${service}-vpce}]" && \
	vpc_endpoint_id=$$(echo $$res | jq -r '.VpcEndpoint.VpcEndpointId') && \
	echo "Saving endpoint id to SSM Parameter Store" && \
	aws ssm put-parameter \
		--name "/${APP_NAME}/${STAGE}/${AWS_REGION}/$${service}EndpointId" \
		--value "$$vpc_endpoint_id" \
		--type String \
		--overwrite
	$(call blue, Successfully created VPC endpoint for $$service)

# test the return code of aws cloudformation deploy and exit if it is 1
service.deploy:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying CloudFormation" 2>&1 | tee -a $${CFN_LOG_PATH}
	@aws cloudformation deploy \
		--region $${AWS_REGION} \
		--template-file template.yaml \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" \
		--tags stage="$${STAGE}" app="$${APP_NAME}" service="${SERVICE}" branch="$$(git branch --show-current)" commit=$$(git rev-parse HEAD) \
		--parameter-overrides \
			Stage="$${STAGE}" \
			AppName="$${APP_NAMEXXX}" \
			createVpc="$${CREATE_VPC}" \
			usePrivateSubnet="$${USE_PRIVATE_SUBNET}" \
			AdminEmail="$${ADMIN_EMAIL}" \
			DataBucketName="$$DATA_BUCKET_NAME" \
			VpcId="$${VPC_ID}" \
			PublicSubnetId="$${PUBLIC_SUBNET_ID}" \
			PrivateSubnetId="$${PRIVATE_SUBNET_ID}" \
			GitHubPersonalAccessToken="$$GITHUB_PERSONAL_ACCESS_TOKEN"

service.monitoring.create-subscriptions:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Creating SNS topic subscriptions" 2>&1 | tee -a $${CFN_LOG_PATH}
	@topic_arn=$$(aws ssm get-parameters \
		--names "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/DataPipelineErrorsTopicArn" \
		--with-decryption \
		--query "Parameters[0].Value" \
		--output text) && \
	for EMAIL in $$(echo $${SUBSCRIBE_EMAILS} | sed 's/,/ /g'); do \
		aws sns subscribe \
			--topic-arn "$$topic_arn" \
			--protocol email \
			--notification-endpoint "$$EMAIL" 2>&1 | tee -a $${CFN_LOG_PATH} || true; \
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
	$(MAKE) service.delete

service.delete.delete-endpoint: #=> service=<string>
	@	if [ -z $$service ]; then \
		echo "No service specified" && exit 1; \
	elif [ $$service != "ssm" ] && [ $$service != "secretsmanager" ]; then \
		echo "Invalid service specified. Valid services are 'ssm', 'secretsmanager'" && exit 1; \
	fi
	@endpoint_id=$$(aws ssm get-parameters \
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