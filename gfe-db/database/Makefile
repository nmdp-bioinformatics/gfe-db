SERVICE := database

target:
	$(info ${HELP_MESSAGE})
	@exit 0

deploy:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying ${SERVICE} service" 2>&1 | tee -a $$CFN_LOG_PATH
	$(MAKE) service.config.deploy
	$(MAKE) service.deploy

service.config.scripts.deploy:
	@config_s3_path=s3://${DATA_BUCKET_NAME}/${CONFIG_S3_PATH}/${SERVICE}/scripts/ && \
	echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying server scripts to $$config_s3_path" 2>&1 \
		| tee -a $$CFN_LOG_PATH && \
	aws s3 cp --recursive --quiet scripts/ $$config_s3_path

service.config.scripts.sync: service.config.scripts.deploy service.config.neo4j.deploy
	@document_name="$$(aws ssm get-parameter \
		--name "/${APP_NAME}/${STAGE}/${AWS_REGION}/DatabaseSyncScriptsDocumentName" | jq -r '.Parameter.Value')" && \
	res=$$(aws ssm send-command \
		--document-name "$$document_name" \
		--targets "Key=instanceids,Values=$${INSTANCE_ID}" \
		--comment "${STAGE}-${APP_NAME} sync scripts to database server" \
		--cloud-watch-output-config '{"CloudWatchOutputEnabled":true}') && \
	command_id=$$(echo "$$res" | jq -r '.Command.CommandId') && \
	echo "\033[0;34mWaiting for command $$command_id to complete on server $${INSTANCE_ID}\033[0m" && \
	aws ssm wait command-executed \
		--command-id "$$command_id" \
		--instance-id "${INSTANCE_ID}" && \
	aws ssm get-command-invocation \
		--command-id "$$command_id" \
		--instance-id "$${INSTANCE_ID}" | jq -r '.StatusDetails' | tee -a $$CFN_LOG_PATH

service.config.cloudwatch-agent.deploy:
	@mkdir -p amazon-cloudwatch-agent/tmp/ 
	@cat amazon-cloudwatch-agent/amazon-cloudwatch-agent.json | sed "s/STAGE/$${STAGE}/g" | sed "s/APP_NAME/$${APP_NAME}/g" > amazon-cloudwatch-agent/tmp/_amazon-cloudwatch-agent.json
	@config_s3_path=s3://${DATA_BUCKET_NAME}/${CONFIG_S3_PATH}/${SERVICE}/amazon-cloudwatch-agent/amazon-cloudwatch-agent.json && \
	echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying amazon-cloudwatch-agent config to $$config_s3_path" 2>&1 \
		| tee -a $$CFN_LOG_PATH && \
	aws s3 cp amazon-cloudwatch-agent/tmp/_amazon-cloudwatch-agent.json $$config_s3_path
	@rm -rf amazon-cloudwatch-agent/tmp/

service.config.neo4j.deploy:
	@aws s3 cp --recursive --quiet neo4j/cypher/ s3://${DATA_BUCKET_NAME}/${CONFIG_S3_PATH}/${SERVICE}/neo4j/cypher/ 2>&1 | tee -a $$CFN_LOG_PATH

service.config.deploy: service.config.scripts.deploy service.config.cloudwatch-agent.deploy service.config.neo4j.deploy

service.ssl.renew-cert:
	@res=$$(aws ssm send-command \
		--document-name "AWS-RunShellScript" \
		--parameters "commands=[cd /home/ec2-user && sudo make ssl.renew-cert]" \
		--targets "Key=instanceids,Values=$${INSTANCE_ID}" \
		--comment "${STAGE}-${APP_NAME} SSL certificate renewal utility" \
		--cloud-watch-output-config '{"CloudWatchOutputEnabled":true}') && \
	command_id=$$(echo "$$res" | jq -r '.Command.CommandId') && \
	aws ssm wait command-executed \
		--command-id "$$command_id" \
		--instance-id "$${INSTANCE_ID}" && \
	aws ssm get-command-invocation \
		--command-id "$$command_id" \
		--instance-id "$${INSTANCE_ID}" | jq -r '.StatusDetails' | tee -a $$CFN_LOG_PATH

service.deploy:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying Neo4j to EC2" 2>&1 | tee -a $$CFN_LOG_PATH
	@aws cloudformation deploy \
		--region $${AWS_REGION} \
		--template-file template.yaml \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" \
		--tags stage="$${STAGE}" app="$${APP_NAME}" service="${SERVICE}" branch="$$(git branch --show-current)" commit=$$(git rev-parse HEAD) \
		--capabilities CAPABILITY_IAM \
		--parameter-overrides \
			Stage="$${STAGE}" \
			AppName="$${APP_NAME}" \
			createVpc="$${CREATE_VPC}" \
			usePrivateSubnet="$${USE_PRIVATE_SUBNET}" \
			DataBucketName="$${DATA_BUCKET_NAME}" \
			Neo4jAmiId="$${NEO4J_AMI_ID}" \
			Neo4jDatabaseName="$${NEO4J_DATABASE_NAME}" \
			Neo4jPassword="$${NEO4J_PASSWORD}" \
			ServiceName="${SERVICE}" \
			ConfigS3Path="$${CONFIG_S3_PATH}" \
			CreateNeo4jUsers="$${CREATE_NEO4J_USERS}" \
			DatabaseVolumeSize="$${DATABASE_VOLUME_SIZE}" \
			HostDomain="$${HOST_DOMAIN}" \
			Subdomain="$${SUBDOMAIN}" \
			AdminEmail="$${ADMIN_EMAIL}" \
			APOCVersion="$${APOC_VERSION}" \
			GDSVersion="$${GDS_VERSION}"

service.backup:
	@res=$$(aws ssm send-command \
		--document-name "AWS-RunShellScript" \
		--parameters "commands=[cd /home/ec2-user && sudo make neo4j.backup]" \
		--targets "Key=instanceids,Values=$${INSTANCE_ID}" \
		--comment "${STAGE}=${APP_NAME} backup service" \
		--cloud-watch-output-config '{"CloudWatchOutputEnabled":true}') && \
	command_id=$$(echo "$$res" | jq -r '.Command.CommandId') && \
	aws ssm wait command-executed \
		--command-id "$$command_id" \
		--instance-id "$${INSTANCE_ID}" && \
	aws ssm get-command-invocation \
		--command-id "$$command_id" \
		--instance-id "$${INSTANCE_ID}" | jq -r '.StatusDetails' | tee -a $$CFN_LOG_PATH

service.backup.list:
	@aws s3 ls --recursive --human-readable s3://${DATA_BUCKET_NAME}/backups/${SERVICE}/${NEO4J_DATABASE_NAME}/

service.backup.get-current:
	@res=$$(aws ssm get-parameter \
		--name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/CurrentBackupS3Path" | jq -r '.Parameter.Value') && \
	echo "$$res"

service.restore: #from_path=s3://<backup path>
	@[ "$$from_path" != "" ] || (echo "from_path is required" && exit 1) && \
	echo "$$from_path" | grep -E '^s3:\/\/[a-zA-Z0-9\.-]+\/([a-zA-Z0-9\.-]+\/)*[a-zA-Z0-9\.-]+\.zip$$' > /dev/null || (echo "\`from_path\` must be an S3 URI: 's3://<bucket>/path'" && exit 1) && \
	echo "Restore $${APP_NAME} from S3 path \`$$from_path\`? [y/N] \c " && read ans && [ $${ans:-N} = y ]
	$(MAKE) service.restore.pre-validate from_path=$$from_path
	@echo "Restoring Neo4j database..." && \
	res=$$(aws ssm send-command \
		--document-name "AWS-RunShellScript" \
		--parameters "commands=[cd /home/ec2-user && sudo make neo4j.restore from_path=$$from_path]" \
		--targets "Key=instanceids,Values=$${INSTANCE_ID}" \
		--comment "${STAGE}=${APP_NAME} restore service" \
		--cloud-watch-output-config '{"CloudWatchOutputEnabled":true}') && \
	command_id=$$(echo "$$res" | jq -r '.Command.CommandId') && \
	aws ssm wait command-executed \
		--command-id "$$command_id" \
		--instance-id "$${INSTANCE_ID}" && \
	aws ssm get-command-invocation \
		--command-id "$$command_id" \
		--instance-id "$${INSTANCE_ID}" | jq -r '.StatusDetails' | tee -a $$CFN_LOG_PATH

service.restore.pre-validate:
	@echo "Validating restore pre-conditions"
	@echo "Checking if backup exists for $$from_path"
	@backup_object=$$(aws s3 ls $$from_path | cut -d ' ' -f 5-) && \
	[ "$$backup_object" != "" ] || (echo "ERROR: Backup does not exist" && exit 1) && \
	echo "Found backup target for $$backup_object"

service.connect:
	@echo "Connecting to public Neo4j EC2 instance"
	@instance_ip=$$(aws ssm get-parameter \
		--name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/Neo4jDatabaseEndpoint" | jq -r '.Parameter.Value') && \
	echo "Connecting to Neo4j EC2 instance at $$instance_ip" && \
	ssh -o "IdentitiesOnly yes" -tt -i ${ROOT_DIR}/${EC2_KEY_PAIR_NAME}.pem ec2-user@$$instance_ip

service.start-session:
	@res=$$(aws ssm get-parameters --names "/${APP_NAME}/${STAGE}/${AWS_REGION}/Neo4jDatabaseInstanceId" | jq -r '.Parameters[0].Value') && \
	aws ssm start-session --target "$$res"

delete: ##=> Delete resources
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting ${SERVICE} service" 2>&1 | tee -a $$CFN_LOG_PATH
	$(MAKE) service.delete

service.delete:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting Neo4j EC2 instance" 2>&1 | tee -a $$CFN_LOG_PATH
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