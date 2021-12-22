SERVICE := infrastructure

target:
	$(info ${HELP_MESSAGE})
	@exit 0

# TODO: Check for secrets in environment variables and cancel if not found
deploy: deploy.cfn

deploy.cfn:
	$(info [*] Deploying ${SERVICE}...)
	@aws cloudformation deploy \
		--region $${REGION} \
		--template-file template.yaml \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" \
		--parameter-overrides \
			Stage="$${STAGE}" \
			AppName="$${APP_NAME}" \
			DataBucketName="$$DATA_BUCKET_NAME" \
			Neo4jUsername="$$NEO4J_USERNAME" \
			Neo4jPassword="$$NEO4J_PASSWORD" \
			GitHubPersonalAccessToken="$$GITHUB_PERSONAL_ACCESS_TOKEN"

delete: ##=> Delete resources
	$(info [*] Deleting infrastructure...)
	@aws cloudformation delete-stack \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" >/dev/null 2>&1 || true && \
	aws cloudformation wait stack-delete-complete \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" >/dev/null 2>&1 || true

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