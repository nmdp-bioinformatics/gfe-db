SERVICE := database
EC2_KEY_PAIR_NAME := $${STAGE}-$${APP_NAME}-$${REGION}-neo4j-key

target:
	$(info ${HELP_MESSAGE})
	@exit 0

deploy:  deploy.check-key-pair deploy.cfn

deploy.check-key-pair: ##=> Checks if the key pair already exists and creates it if it does not
	@key_pair="$$(aws ec2 describe-key-pairs --key-name ${EC2_KEY_PAIR_NAME} | jq '.KeyPairs[0].KeyName' || true)" && \
	[ "$$key_pair" ] && echo "Key pair found: $$key_pair" && exit 0 || echo "No key pair found" && \
	echo "Creating EC2 key pair \"${EC2_KEY_PAIR_NAME}\"..." && \
	aws ec2 create-key-pair --key-name ${EC2_KEY_PAIR_NAME} | jq -r '.KeyMaterial' > ${EC2_KEY_PAIR_NAME}.pem && \
	aws ssm put-parameter \
		--name "/$${APP_NAME}/$${STAGE}/$${REGION}/EC2KeyPairName" \
		--type "String" \
		--value "${EC2_KEY_PAIR_NAME}" \
		--overwrite >/dev/null 2>&1

deploy.cfn:
	$(info [*] Deploying ${SERVICE}...)
	@aws cloudformation deploy \
		--region $${REGION} \
		--template-file template.yaml \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" \
		--capabilities CAPABILITY_IAM \
		--parameter-overrides \
			Stage="$${STAGE}" \
			AppName="$${APP_NAME}"

delete: ##=> Delete resources
	$(info [*] Deleting database...)
	@aws cloudformation delete-stack \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" >/dev/null 2>&1 || true && \
	aws cloudformation wait stack-delete-complete \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" >/dev/null 2>&1 || true && \
	aws ec2 delete-key-pair --key-name "${EC2_KEY_PAIR_NAME}" >/dev/null 2>&1 || true && \
	aws ssm delete-parameter --name "/$${APP_NAME}/$${STAGE}/$${REGION}/EC2KeyPairName" >/dev/null 2>&1 || true

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