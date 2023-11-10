SERVICE := access

target:
	$(info ${HELP_MESSAGE})
	@exit 0

# TODO parameterize IP address for security group
deploy:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying ${SERVICE} service" 2>&1 | tee -a $$CFN_LOG_PATH
	$(MAKE) service.deploy

service.deploy: 
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying Neo4j to EC2 with ${NEO4J_AMI_ID}" 2>&1 | tee -a $$CFN_LOG_PATH
	@aws cloudformation deploy \
		--region $${AWS_REGION} \
		--template-file template.yaml \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" \
		--tags stage="$${STAGE}" app="$${APP_NAME}" service="${SERVICE}" branch="$$(git branch --show-current)" commit=$$(git rev-parse HEAD) \
		--capabilities CAPABILITY_IAM \
		--parameter-overrides \
			Stage="$${STAGE}" \
			AppName="$${APP_NAME}" \
			Ec2KeyPairName="${EC2_KEY_PAIR_NAME}" \
			MyIp="0.0.0.0/0" 2>&1
	$(MAKE) service.authorize-security-group-ingress

# connect to the bastion host using ssh
# 1) get the bastion host's public ip address from SSM Parameter Store using the path /${AppName}/${Stage}/${AWS::Region}/BastionHostPublicIp
# 2) connect to the bastion host using ssh
# service.connect:
# 	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Connecting to Neo4j database using bastion host" 2>&1 | tee -a $$CFN_LOG_PATH
# 	@bastion_host_public_ip="$$(aws ssm get-parameter \
# 		--name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/BastionHostPublicIp" | jq -r '.Parameter.Value')" && \
# 	ssh -i ${ROOT_DIR}/${EC2_KEY_PAIR_NAME}.pem ec2-user@"$$bastion_host_public_ip"

service.bastion-server.connect:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Setting up SSH tunnel to database instance using bastion host" 2>&1 | tee -a $$CFN_LOG_PATH
	@bastion_host_public_ip="$$(aws ssm get-parameter \
		--name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/BastionHostPublicIp" | jq -r '.Parameter.Value')" && \
	database_private_ip="$$(aws ssm get-parameter \
		--name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/Neo4jPrivateIp" | jq -r '.Parameter.Value')" && \
	ssh-add ${ROOT_DIR}/${EC2_KEY_PAIR_NAME}.pem && \
	ssh -tt -i ${ROOT_DIR}/${EC2_KEY_PAIR_NAME}.pem -A ec2-user@"$$bastion_host_public_ip" ssh bitnami@"$$database_private_ip"


# TODO add the bastion host's security group to the neo4j security group
# 1) get the bastion host's security group id from SSM Parameter Store using the path /${AppName}/${Stage}/${AWS::Region}/BastionHostSecurityGroupId
# 2) get the neo4j security group id from SSM Parameter Store using the path /${AppName}/${Stage}/${AWS::Region}/Neo4jDatabaseSecurityGroupId
# 3) add the bastion host's security group id to the neo4j security group
service.authorize-security-group-ingress:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Updating Neo4j database security group for bastion host" 2>&1 | tee -a $$CFN_LOG_PATH
	@bastion_host_security_group_id="$$(aws ssm get-parameter \
		--name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/BastionHostSecurityGroupId" | jq -r '.Parameter.Value')" && \
	neo4j_database_security_group_id="$$(aws ssm get-parameter \
		--name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/Neo4jDatabaseSecurityGroupId" | jq -r '.Parameter.Value')" && \
	aws ec2 authorize-security-group-ingress \
		--group-id "$$neo4j_database_security_group_id" \
		--protocol tcp \
		--port 22 \
		--source-group "$$bastion_host_security_group_id" 2>&1 | tee -a $$CFN_LOG_PATH || true

delete: ##=> Delete resources
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting ${SERVICE} service" 2>&1 | tee -a $$CFN_LOG_PATH
	$(MAKE) service.delete

service.revoke-security-group-ingress:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Removing bastion host security group from Neo4j database security group" 2>&1 | tee -a $$CFN_LOG_PATH
	@bastion_host_security_group_id="$$(aws ssm get-parameter \
		--name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/BastionHostSecurityGroupId" | jq -r '.Parameter.Value')" && \
	neo4j_database_security_group_id="$$(aws ssm get-parameter \
		--name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/Neo4jDatabaseSecurityGroupId" | jq -r '.Parameter.Value')" && \
	aws ec2 revoke-security-group-ingress \
		--group-id "$$neo4j_database_security_group_id" \
		--protocol tcp \
		--port 22 \
		--source-group "$$bastion_host_security_group_id" 2>&1 | tee -a $$CFN_LOG_PATH || true

service.delete: service.revoke-security-group-ingress
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting Bastion EC2 instance" 2>&1 | tee -a $$CFN_LOG_PATH
	@aws cloudformation delete-stack \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" 2>&1 | tee -a $$CFN_LOG_PATH || true && \
	aws cloudformation wait stack-delete-complete \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" 2>&1 | tee -a $$CFN_LOG_PATH || true

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