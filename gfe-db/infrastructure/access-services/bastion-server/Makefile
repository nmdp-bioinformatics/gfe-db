SERVICE := bastion-server

target:
	$(info ${HELP_MESSAGE})
	@exit 0

# TODO parameterize IP address for security group
deploy:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying ${SERVICE} service" 2>&1 | tee -a $$CFN_LOG_PATH
	$(MAKE) service.deploy

service.deploy: #=> admin_ip=<ip address>
	@echo "\033[0;34m$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying Bastion Server\033[0m" 2>&1 | tee -a $${CFN_LOG_PATH}
	@echo "\033[0;32mFound Admin IP: $${ADMIN_IP}\033[0m"
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
			AdminIp="$${ADMIN_IP}"
	$(MAKE) service.authorize-security-group-ingress

service.connect:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Setting up SSH tunnel to database instance using bastion host" 2>&1 | tee -a $$CFN_LOG_PATH
	@bastion_host_public_ip="$$(aws ssm get-parameter \
		--name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/BastionHostElasticIp" | jq -r '.Parameter.Value')" && \
	echo "\033[0;32mBastion host public IP: $$bastion_host_public_ip\033[0m" && \
	database_private_ip="$$(aws ssm get-parameter \
		--name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/Neo4jPrivateIp" | jq -r '.Parameter.Value')" && \
	echo "\033[0;32mDatabase private IP: $$database_private_ip\033[0m" && \
	ssh-add -K ${ROOT_DIR}/${EC2_KEY_PAIR_NAME}.pem && \
	ssh -tt -o IdentitiesOnly=yes -i ${ROOT_DIR}/${EC2_KEY_PAIR_NAME}.pem -A ec2-user@"$$bastion_host_public_ip" ssh ubuntu@"$$database_private_ip"

service.ui.connect:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Setting up SSH tunnel to database instance using bastion host" 2>&1 | tee -a $$CFN_LOG_PATH
	@bastion_host_public_ip="$$(aws ssm get-parameter \
		--name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/BastionHostElasticIp" | jq -r '.Parameter.Value')" && \
	database_private_ip="$$(aws ssm get-parameter \
		--name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/Neo4jPrivateIp" | jq -r '.Parameter.Value')" && \
	echo "\033[0;32mNeo4j Browser is available at: http://localhost:7474/browser/\033[0m" && \
	ssh -L 7474:$$database_private_ip:7474 -L 7687:$$database_private_ip:7687 -N -J ec2-user@$$bastion_host_public_ip ubuntu@$$database_private_ip

service.authorize-security-group-ingress:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Updating Neo4j database security group for bastion host" 2>&1 | tee -a $$CFN_LOG_PATH
	@bastion_host_security_group_id="$$(aws ssm get-parameter \
			--name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/BastionHostSecurityGroupId" | jq -r '.Parameter.Value')" && \
	echo "\033[0;32mFound bastion host security group ID: $$bastion_host_security_group_id\033[0m" && \
	neo4j_database_security_group_id="$$(aws ssm get-parameter \
		--name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/Neo4jDatabaseSecurityGroupId" | jq -r '.Parameter.Value')" && \
	echo "\033[0;32mFound Neo4j database security group ID: $$neo4j_database_security_group_id\033[0m" && \
	echo "\033[0;34mValidating security group configuration...\033[0m" && \
	res=$$(aws ec2 describe-security-groups \
		--group-ids "$$neo4j_database_security_group_id" \
		| jq -r '.SecurityGroups[].IpPermissions[] | select(.FromPort == 22 and .ToPort == 22 and .IpProtocol == "tcp") | .UserIdGroupPairs[].GroupId') && \
	if [ "$$res" != "$$bastion_host_security_group_id" ]; then \
		echo "\033[0;34mAuthorizing bastion host access to database...\033[0m"; \
		res=$$(aws ec2 authorize-security-group-ingress \
			--group-id "$$neo4j_database_security_group_id" \
			--protocol tcp \
			--port 22 \
			--source-group "$$bastion_host_security_group_id") && \
		echo "\033[0;34mBastion host security group is now associated with Neo4j database security group\033[0m" || \
		echo "\033[0;31mFailed to associate bastion host security group with Neo4j database security group\033[0m"; \
	else \
		echo "\033[0;34mBastion host security group is already associated with Neo4j database security group\033[0m"; \
	fi

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

	...::: Connect to database instance using bastion host :::...
	$ make service.connect

endef