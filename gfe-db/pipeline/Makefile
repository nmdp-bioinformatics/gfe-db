SERVICE := pipeline

target:
	$(info ${HELP_MESSAGE})
	@exit 0

deploy: deploy.cfn deploy.jobs deploy.functions

deploy.cfn:
	$(info [*] Deploying ${SERVICE}...)
	@aws cloudformation deploy \
		--region $${REGION} \
		--template-file template.yaml \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" \
		--capabilities CAPABILITY_IAM \
		--parameter-overrides \
			Stage="$${STAGE}" \
			AppName="$${APP_NAME}" \
			BuildServiceRepositoryName="${BUILD_REPOSITORY}" \
			LoadServiceRepositoryName="${LOAD_REPOSITORY}"

deploy.jobs:
	$(MAKE) -C jobs/ deploy

deploy.functions:
	$(MAKE) -C functions/ deploy

delete: delete.functions delete.jobs delete.cfn

delete.cfn:
	$(info [*] Deleting resources...)
	@aws cloudformation delete-stack \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" >/dev/null 2>&1

delete.jobs:
	$(MAKE) -C jobs/ delete

delete.functions:
	$(MAKE) -C functions/ delete


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