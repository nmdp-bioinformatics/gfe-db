SERVICE := s3-backup

target:
	$(info ${HELP_MESSAGE})
	@exit 0

deploy: deploy.backup

deploy.backup:
	$(info [*] Deploying ${SERVICE}...)
	@aws cloudformation deploy \
		--region $${REGION} \
		--template-file template.yaml \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" \
		--capabilities CAPABILITY_IAM \
		--parameter-overrides \
			Stage="$${STAGE}" \
			AppName="$${APP_NAME}"
	@echo "Deploying backup script to S3..." && \
	aws s3 cp backup.sh s3://$$DATA_BUCKET_NAME/config/scripts/backup.sh

delete: delete.backup ##=> Delete resources

delete.backup:
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