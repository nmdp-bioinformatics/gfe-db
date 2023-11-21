SERVICE := database

target:
	$(info ${HELP_MESSAGE})
	@exit 0

deploy:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying ${SERVICE} service" 2>&1 | tee -a $$CFN_LOG_PATH
	$(MAKE) service.config.scripts.deploy
	$(MAKE) service.config.cloudwatch-agent.deploy
	$(MAKE) service.deploy
	$(MAKE) service.config.neo4j.deploy

service.config.neo4j.deploy:
ifeq ($(USE_PRIVATE_SUBNET),true)
	@aws s3 cp neo4j/neo4j.private.template s3://$$DATA_BUCKET_NAME/config/neo4j/neo4j.conf 2>&1 | tee -a $$CFN_LOG_PATH
else ifeq ($(USE_PRIVATE_SUBNET),false)
	@aws s3 cp neo4j/neo4j.template s3://$$DATA_BUCKET_NAME/config/neo4j/neo4j.conf 2>&1 | tee -a $$CFN_LOG_PATH
endif
	@aws s3 cp --recursive neo4j/cypher/ s3://$$DATA_BUCKET_NAME/config/neo4j/cypher/ 2>&1 | tee -a $$CFN_LOG_PATH

service.config.scripts.deploy:
	@script_s3_path=s3://$$DATA_BUCKET_NAME/config/scripts/ && \
	echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying backup script to $$script_s3_path" 2>&1 \
		| tee -a $$CFN_LOG_PATH && \
	aws s3 cp --recursive scripts/ $$script_s3_path

service.config.scripts.sync: service.config.scripts.deploy
	@document_name="$$(aws ssm get-parameter \
		--name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/DatabaseSyncScriptsDocumentName" | jq -r '.Parameter.Value')" && \
	aws ssm send-command \
		--document-name "$$document_name" \
		--targets "Key=instanceids,Values=${INSTANCE_ID}" \
		--comment "${STAGE}=${APP_NAME} backup service" 2>&1 | tee -a $$CFN_LOG_PATH || true

service.ssl.renew-cert:
	@aws ssm send-command \
		--document-name "AWS-RunShellScript" \
		--parameters "commands=[cd /home/bitnami && sudo make ssl.renew-cert]" \
		--targets "Key=instanceids,Values=${INSTANCE_ID}" \
		--comment "${STAGE}-${APP_NAME} SSL certificate renewal utility" 2>&1 | tee -a $$CFN_LOG_PATH || true

service.config.cloudwatch-agent.deploy:
	@mkdir -p amazon-cloudwatch-agent/tmp/ 
	@cat amazon-cloudwatch-agent/amazon-cloudwatch-agent.json | sed "s/STAGE/$${STAGE}/g" | sed "s/APP_NAME/$${APP_NAME}/g" > amazon-cloudwatch-agent/tmp/_amazon-cloudwatch-agent.json
	@config_s3_path=s3://$$DATA_BUCKET_NAME/config/amazon-cloudwatch-agent/amazon-cloudwatch-agent.json && \
	echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying amazon-cloudwatch-agent config to $$config_s3_path" 2>&1 \
		| tee -a $$CFN_LOG_PATH && \
	aws s3 cp amazon-cloudwatch-agent/tmp/_amazon-cloudwatch-agent.json $$config_s3_path
	@rm -rf amazon-cloudwatch-agent/tmp/

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
			createVpc="$${CREATE_VPC}" \
			usePrivateSubnet="$${USE_PRIVATE_SUBNET}" \
			DataBucketName="$${DATA_BUCKET_NAME}" \
			Neo4jAmiId="$${NEO4J_AMI_ID}" \
			Neo4jPassword="$${NEO4J_PASSWORD}" \
			DatabaseVolumeSize="$${DATABASE_VOLUME_SIZE}" \
			HostDomain="$${HOST_DOMAIN}" \
			Subdomain="$${SUBDOMAIN}" \
			AdminEmail="$${ADMIN_EMAIL}" \
			APOCVersion="$${APOC_VERSION}" \
			GDSVersion="$${GDS_VERSION}" 2>&1 | tee -a $$CFN_LOG_PATH || true

# TODO repurpose for updating neo4j.conf
service.config.neo4j.update: service.config.update-dns.private service.config.neo4j.deploy
	@aws ssm send-command \
		--document-name "AWS-RunShellScript" \
		--parameters "commands=[cd /home/bitnami && sudo make neo4j.config.update]" \
		--targets "Key=instanceids,Values=${INSTANCE_ID}" \
		--comment "${STAGE}=${APP_NAME} backup service" 2>&1 | tee -a $$CFN_LOG_PATH || true

service.backup:
	@aws ssm send-command \
		--document-name "AWS-RunShellScript" \
		--parameters "commands=[cd /home/bitnami && sudo make neo4j.backup]" \
		--targets "Key=instanceids,Values=${INSTANCE_ID}" \
		--comment "${STAGE}=${APP_NAME} backup service" 2>&1 | tee -a $$CFN_LOG_PATH || true

service.backup.list:
	@aws s3 ls --recursive s3://${DATA_BUCKET_NAME}/backups/neo4j/ | cut -d' ' -f7 | cut -d'/' -f3-6

service.restore: #from_path=s3://<backup path>
	@[ "$$from_path" != "" ] || (echo "from_path is required" && exit 1) && \
	echo "$$from_path" | grep -E '^s3:\/\/[a-zA-Z0-9\.-]+\/([a-zA-Z0-9\.-]+\/)*[a-zA-Z0-9\.-]+\.zip$$' > /dev/null || (echo "\`from_path\` must be an S3 URI: 's3://<bucket>/path'" && exit 1) && \
	echo "Restore $${APP_NAME} from S3 path \`$$from_path\`? [y/N] \c " && read ans && [ $${ans:-N} = y ]
	$(MAKE) service.restore.pre-validate from_path=$$from_path
	@aws ssm send-command \
		--document-name "AWS-RunShellScript" \
		--parameters "commands=[cd /home/bitnami && sudo make neo4j.restore from_path=$$from_path]" \
		--targets "Key=instanceids,Values=${INSTANCE_ID}" \
		--comment "${STAGE}=${APP_NAME} restore service" 2>&1 | tee -a $$CFN_LOG_PATH || true

service.restore.pre-validate:
	@echo "Validating restore pre-conditions"
	@echo "Checking if backup exists for $$from_path"
	@backup_object=$$(aws s3 ls $$from_path | cut -d ' ' -f 6-) && \
	[ "$$backup_object" != "" ] || (echo "ERROR: Backup does not exist" && exit 1) && \
	echo "Found backup target for $$backup_object"

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