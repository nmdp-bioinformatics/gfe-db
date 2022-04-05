SERVICE := database
EC2_KEY_PAIR_NAME := $${STAGE}-$${APP_NAME}-$${REGION}-neo4j-key

target:
	$(info ${HELP_MESSAGE})
	@exit 0

deploy:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying ${SERVICE} service" 2>&1 | tee -a $$CFN_LOG_PATH
	$(MAKE) deploy.check-key-pair
	$(MAKE) deploy.config
	$(MAKE) deploy.backup
	$(MAKE) deploy.database
	$(MAKE) deploy.copy-s3-logs


deploy.check-key-pair: ##=> Checks if the key pair already exists and creates it if it does not
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Checking for key pair ${EC2_KEY_PAIR_NAME}" 2>&1 | tee -a $$CFN_LOG_PATH
	@key_pair="$$(aws ec2 describe-key-pairs --key-name ${EC2_KEY_PAIR_NAME} | jq '.KeyPairs[0].KeyName' || true)" && \
	[ "$$key_pair" ] && echo "Key pair found: $$key_pair" && exit 0 || echo "No key pair found..." && \
	echo "Creating EC2 key pair \"${EC2_KEY_PAIR_NAME}\"" && \
	aws ec2 create-key-pair --key-name ${EC2_KEY_PAIR_NAME} | jq -r '.KeyMaterial' > ${EC2_KEY_PAIR_NAME}.pem && \
	aws ssm put-parameter \
		--name "/$${APP_NAME}/$${STAGE}/$${REGION}/EC2KeyPairName" \
		--type "String" \
		--value "${EC2_KEY_PAIR_NAME}" \
		--overwrite 2>&1 | tee -a $$CFN_LOG_PATH || true

deploy.config: config.neo4j

config.neo4j:
	@neo4j_conf_s3_path=s3://$$DATA_BUCKET_NAME/config/neo4j/neo4j.template && \
	echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying Neo4j config to $$neo4j_conf_s3_path" 2>&1 \
		| tee -a $$CFN_LOG_PATH && \
	aws s3 cp neo4j/neo4j.template $$neo4j_conf_s3_path 2>&1 | tee -a $$CFN_LOG_PATH
	@neo4j_cypher_s3_path=s3://$$DATA_BUCKET_NAME/config/neo4j/cypher/ && \
	echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying Neo4j config to $$neo4j_cypher_s3_path" 2>&1 \
		| tee -a $$CFN_LOG_PATH && \
	aws s3 cp --recursive neo4j/cypher/ $$neo4j_cypher_s3_path 2>&1 | tee -a $$CFN_LOG_PATH

deploy.backup:
	@script_s3_path=s3://$$DATA_BUCKET_NAME/config/scripts/backup.sh && \
	echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying backup script to $$script_s3_path" 2>&1 \
		| tee -a $$CFN_LOG_PATH && \
	aws s3 cp scripts/backup.sh $$script_s3_path

# TODO: get endpoint and store as SSM parameter to replace resource in CloudFormation template
deploy.database:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying Neo4j to EC2 with ${NEO4J_AMI_ID}" 2>&1 | tee -a $$CFN_LOG_PATH
	@aws cloudformation deploy \
		--region $${REGION} \
		--template-file template.yaml \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" \
		--capabilities CAPABILITY_IAM \
		--parameter-overrides \
			Stage="$${STAGE}" \
			AppName="$${APP_NAME}" \
			DataBucketName="$$DATA_BUCKET_NAME" \
			Neo4jAmiId="${NEO4J_AMI_ID}" \
			Neo4jPassword="$$NEO4J_PASSWORD" \
			DatabaseVolumeSize="$$DATABASE_VOLUME_SIZE" 2>&1 | tee -a $$CFN_LOG_PATH || true


deploy.copy-s3-logs:
	@aws s3 cp --recursive s3://$$DATA_BUCKET_NAME/logs/ $$LOGS_DIR/

# ec2.check-state:
# 	@instance_id="$$(aws ssm get-parameters \
# 		--names "/$${APP_NAME}/$${STAGE}/$${REGION}/Neo4jDatabaseInstanceId" \
# 			| jq -r '.Parameters | map(select(.Version == 1))[0].Value')" && \
# 	echo $$instance_id && \
# 	instance_status="$$(aws ec2 describe-instance-status --instance-ids $$instance_id)" | jq && \
# 	echo $$instance_status

delete: ##=> Delete resources
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting ${SERVICE} service" 2>&1 | tee -a $$CFN_LOG_PATH
	$(MAKE) delete.database

delete.database:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting Neo4j EC2 instance" 2>&1 | tee -a $$CFN_LOG_PATH
	@aws cloudformation delete-stack \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" 2>&1 | tee -a $$CFN_LOG_PATH || true && \
	aws cloudformation wait stack-delete-complete \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" 2>&1 | tee -a $$CFN_LOG_PATH || true && \
	aws ec2 delete-key-pair --key-name "${EC2_KEY_PAIR_NAME}" 2>&1 | tee -a $$CFN_LOG_PATH || true && \
	aws ssm delete-parameter --name "/$${APP_NAME}/$${STAGE}/$${REGION}/EC2KeyPairName" 2>&1 | tee -a $$CFN_LOG_PATH || true

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