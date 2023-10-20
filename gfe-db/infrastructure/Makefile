.PHONY: service.update-vpc-endpoints
SERVICE := infrastructure

target:
	$(info ${HELP_MESSAGE})
	@exit 0

# service.env.validate:

deploy: 
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying ${SERVICE} service" 2>&1 | tee -a $${CFN_LOG_PATH}
	$(MAKE) service.deploy
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
	echo "Response:" && \
	echo $$res | jq -r

# TODO add cloudformation waiter for stack create complete and exit if failed
# TODO test AWSCLI output and validate the stack was created successfully
# Add stateful check from SSM Param of vpc=true/false to correctly set CREATE_VPC
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
			createSsmVpcEndpoint="$${CREATE_SSM_ENDPOINT}" \
			createSecretsManagerVpcEndpoint="$${CREATE_SECRETSMANAGER_ENDPOINT}" \
			AdminEmail="$${ADMIN_EMAIL}" \
			DataBucketName="$$DATA_BUCKET_NAME" \
			GitHubPersonalAccessToken="$$GITHUB_PERSONAL_ACCESS_TOKEN" \
			HostedZoneId="$${HOSTED_ZONE_ID}" \
			VpcId="$${VPC_ID}" \
			PublicSubnetId="$${PUBLIC_SUBNET_ID}" \
			PrivateSubnetId="$${PRIVATE_SUBNET_ID}"
	$(MAKE) service.update-vpc-endpoints

# Not used for private subnet deployment
# $(MAKE) service.deploy.update-dns

service.update-vpc-endpoints:
ifeq ($(CREATE_SSM_ENDPOINT),false)
	$(MAKE) service.update-vpc-endpoint service=ssm
endif
ifeq ($(CREATE_SECRETSMANAGER_ENDPOINT),false)
	$(MAKE) service.update-vpc-endpoint service=secretsmanager
endif

# TODO get endpoint id from env vars, --filter "Name=subnet-id,Values=$${PRIVATE_SUBNET_ID}" does not work
# the vpc endpoint must match PRIVATE_SUBNET_ID env var
service.update-vpc-endpoint: #service=<string[ssm,secretsmanager]>
	@echo "\033[0;34m$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Adding security group to $${service} VPC endpoint\033[0m" 2>&1 | tee -a $${CFN_LOG_PATH}
	@endpoint_json=$$(aws ec2 describe-vpc-endpoints | jq 'del(.VpcEndpoints[].PolicyDocument)') && \
	endpoint_id=$$(echo $$endpoint_json | jq -r --arg subnet "$${PRIVATE_SUBNET_ID}" '.VpcEndpoints[] | select(.SubnetIds[] == $$subnet) .VpcEndpointId') && \
	if [ -z "$$endpoint_id" ]; then \
		echo "\033[0;31m$${service} VPC endpoint is not associated with \`PRIVATE_SUBNET_ID\` ($${PRIVATE_SUBNET_ID})\033[0m" && exit 1; \
	fi; \
	[ "$${service}" = "ssm" ] && endpoint_id=$$SSM_ENDPOINT_ID; \
	[ "$${service}" = "secretsmanager" ] && endpoint_id=$$SECRETSMANAGER_ENDPOINT_ID; \
	echo "\033[0;32mFound VPC endpoint ID: $$endpoint_id" && \
	security_group_id=$$(aws ssm get-parameters \
		--names "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/Neo4jDatabaseSecurityGroupId" \
		--with-decryption \
		--query "Parameters[0].Value" \
		--output text) && \
	echo "Found security group ID: $$security_group_id" && \
	response=$$(aws ec2 modify-vpc-endpoint \
		--vpc-endpoint-id "$$endpoint_id" \
		--add-security-group-ids "$$security_group_id" | jq -r '.Return') && \
	[ "$$response" = "true" ] && echo "\033[0;34mSuccessfully added security group to $${service} VPC endpoint\033[0m" || \
	(echo "\033[0;31mFailed to add security group to $${service} VPC endpoint\033[0m" && exit 1)

service.remove-security-groups:
ifeq ($(CREATE_SSM_ENDPOINT),false)
	$(MAKE) service.remove-security-group service=ssm
endif
ifeq ($(CREATE_SECRETSMANAGER_ENDPOINT),false)
	$(MAKE) service.remove-security-group service=secretsmanager
endif

# Removes the Neo4j security group from the VPC endpoints
# aws ec2 modify-vpc-endpoint --remove-security-group-ids 
service.remove-security-group: #service=<string[ssm,secretsmanager]>
	@echo "\033[0;34m$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Removing security group from $${service} VPC endpoint\033[0m" 2>&1 | tee -a $${CFN_LOG_PATH}
	@[ "$${service}" = "ssm" ] && endpoint_id=$$SSM_ENDPOINT_ID; \
	[ "$${service}" = "secretsmanager" ] && endpoint_id=$$SECRETSMANAGER_ENDPOINT_ID; \
	echo "\033[0;32mFound VPC endpoint ID: $$endpoint_id\033[0m" && \
	security_group_id=$$(aws ssm get-parameters \
		--names "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/Neo4jDatabaseSecurityGroupId" \
		--with-decryption \
		--query "Parameters[0].Value" \
		--output text) && \
	echo "\033[0;32mFound security group ID: $$security_group_id\033[0m" && \
	response=$$(aws ec2 modify-vpc-endpoint \
		--vpc-endpoint-id "$$endpoint_id" \
		--remove-security-group-ids "$$security_group_id" | jq -r '.Return') && \
	[ "$$response" = "true" ] && echo "\033[0;34mSuccessfully removed security group from $${service} VPC endpoint\033[0m" || \
	(echo "\033[0;31mFailed to remove security group from $${service} VPC endpoint\033[0m" && exit 1) | jq -r

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

# TODO remove the sg association from the vpc endpoints before deleting
delete: ##=> Delete resources
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting ${SERVICE} service" 2>&1 | tee -a $${CFN_LOG_PATH}
ifeq ($(data),true)
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting data in $$DATA_BUCKET_NAME" 2>&1 | tee -a $${CFN_LOG_PATH}
	@aws s3 rm --recursive s3://$${DATA_BUCKET_NAME}/
endif
	$(MAKE) service.remove-security-groups
	$(MAKE) service.delete

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