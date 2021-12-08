SERVICE := infrastructure

target:
	$(info ${HELP_MESSAGE})
	@exit 0

deploy:  deploy.cfn

deploy.cfn:
	$(info [*] Deploying ${SERVICE}...)
	@aws cloudformation deploy \
		--region $${REGION} \
		--template-file template.yaml \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" \
		--parameter-overrides \
			Stage="$${STAGE}" \
			AppName="$${APP_NAME}" \
			EC2KeyPairName="dev-gfe-db-us-east-1-ec2-key" \
			Neo4jUsername="$$NEO4J_USERNAME" \
			Neo4jPassword="$$NEO4J_PASSWORD" \
			GitHubPersonalAccessToken="$$GITHUB_PERSONAL_ACCESS_TOKEN"

delete: ##=> Delete resources
	$(info [*] Deleting resources...)
	@aws cloudformation delete-stack \
		--stack-name $${STACK_NAME}-"$$(cat environment.json | jq -r '."$(environment)".Name')"
	$(info [*] Deleting script in S3...)
	@aws s3 rm --quiet \
		"$$(cat environment.json | jq -r '."$(environment)"."Command.ScriptLocation"')"

#############
#  Helpers  #
#############

define HELP_MESSAGE

	Environment variables:

	JOB_NAME: "$$(cat environment.json | jq -r '."$(environment)".Name')"
		Description: Stack Name already deployed
	REGION: "${REGION}"
		Description: AWS region

	Common usage:

	...::: Deploy all CloudFormation based services :::...
	$ make deploy

endef