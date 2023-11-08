SERVICE := database
export EC2_KEY_PAIR_NAME := $${STAGE}-$${APP_NAME}-$${AWS_REGION}-neo4j-key

target:
	$(info ${HELP_MESSAGE})
	@exit 0

deploy:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying ${SERVICE} service" 2>&1 | tee -a $$CFN_LOG_PATH
	$(MAKE) service.check-key-pair
	$(MAKE) service.config.scripts.deploy
	$(MAKE) service.config.cloudwatch-agent.deploy
	$(MAKE) service.deploy
	$(MAKE) service.config.neo4j.deploy

service.check-key-pair: ##=> Checks if the key pair already exists and creates it if it does not
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Checking for key pair ${EC2_KEY_PAIR_NAME}" 2>&1 | tee -a $$CFN_LOG_PATH
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

service.config.neo4j.deploy:
ifeq ($(USE_PRIVATE_SUBNET),true)
	@aws s3 cp neo4j/neo4j.private.template s3://$$DATA_BUCKET_NAME/config/neo4j/neo4j.conf 2>&1 | tee -a $$CFN_LOG_PATH
else ifeq ($(USE_PRIVATE_SUBNET),false)
	@aws s3 cp neo4j/neo4j.template s3://$$DATA_BUCKET_NAME/config/neo4j/neo4j.conf 2>&1 | tee -a $$CFN_LOG_PATH
endif
	@aws s3 cp --recursive neo4j/cypher/ s3://$$DATA_BUCKET_NAME/config/neo4j/cypher/ 2>&1 | tee -a $$CFN_LOG_PATH

# # This target will configure the database to use the correct DNS name using env vars
# service.config.update-dns.public:
# 	@[ "${HOST_DOMAIN}" != "" ] && \
# 	cat neo4j/neo4j.template | \
# 	sed s"/# dbms.default_advertised_address=.*/dbms.default_advertised_address=${SUBDOMAIN}.${HOST_DOMAIN}/" > neo4j/neo4j.conf

# # TODO Neo4jPrivateIp isn't available when this is run, so it needs to be refactored
# # This target will configure the database to use the correct DNS name using env vars
# service.config.update-dns.private:
# 	@private_ip="$$(aws ssm get-parameter \
# 		--name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/Neo4jPrivateIp" | jq -r '.Parameter.Value')" && \
# 	echo "Private IP: $$private_ip" && \
# 	cat neo4j/neo4j.private.template | \
# 	sed "s|# dbms.default_advertised_address=.*|dbms.default_advertised_address="$$private_ip"|" > neo4j/neo4j.conf

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

service.bastion-server.deploy:
	$(MAKE) -C bastion deploy

service.bastion-server.connect:
	$(MAKE) -C bastion service.connect

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

# TODO add confirmation with prequisites (warnings to avoid losing data, make sure scripts are synced etc.) ✅
# TODO database restore: validate date parameter is given ✅ and in the correct format ✅
# TODO move validation to the database/Makefile ✅
# TODO database restore: validate that the backup exists
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
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" 2>&1 | tee -a $$CFN_LOG_PATH || true && \
	aws ec2 delete-key-pair --key-name "${EC2_KEY_PAIR_NAME}" 2>&1 | tee -a $$CFN_LOG_PATH || true && \
	aws ssm delete-parameter --name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/EC2KeyPairName" 2>&1 | tee -a $$CFN_LOG_PATH || true && \
	mv ${ROOT_DIR}/${EC2_KEY_PAIR_NAME}.pem ${ROOT_DIR}/deprecated-key-$${INSTANCE_ID}-$$(gdate -u +'%Y-%m-%d-%H-%M').pem || true && \
	echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting key pair ${EC2_KEY_PAIR_NAME} in Secrets Manager" 2>&1 | tee -a $$CFN_LOG_PATH && \
	aws secretsmanager delete-secret \
		--secret-id "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/EC2KeyPair" \
		--force-delete-without-recovery 2>&1 | tee -a $$CFN_LOG_PATH || true


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