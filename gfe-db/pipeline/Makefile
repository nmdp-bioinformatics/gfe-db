SERVICE := pipeline

target:
	$(info ${HELP_MESSAGE})
	@exit 0

# TODO: Don't deploy jobs if pipeline stack fails to create (exit Make)
deploy:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying ${SERVICE} service" 2>&1 | tee -a $$CFN_LOG_PATH
	$(MAKE) service.functions.deploy
	$(MAKE) service.jobs.deploy
	$(MAKE) service.config.deploy

service.functions.deploy:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying ${SERVICE} - functions" 2>&1 | tee -a $$CFN_LOG_PATH
	@sam build \
		--template-file template.yaml  2>&1 | tee -a $$CFN_LOG_PATH || true && \
	sam package \
		--resolve-s3 \
		--output-template-file packaged.yaml  2>&1 | tee -a $$CFN_LOG_PATH || true && \
	sam deploy \
		--no-fail-on-empty-changeset \
		--region "$${AWS_REGION}" \
		--template-file packaged.yaml \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" \
		--capabilities CAPABILITY_IAM \
		--parameter-overrides \
			Stage="$${STAGE}" \
			AppName="$${APP_NAME}" \
			InvokePipelineFunctionName="$${STAGE}"-"$${APP_NAME}"-"$$(cat functions/environment.json | jq -r '.Functions.InvokePipeline.FunctionConfiguration.FunctionName')" \
			InvokePipelineFunctionMemorySize="$$(cat functions/environment.json | jq -r '.Functions.InvokePipeline.FunctionConfiguration.MemorySize')" \
			InvokePipelineFunctionTimeout="$$(cat functions/environment.json | jq -r '.Functions.InvokePipeline.FunctionConfiguration.Timeout')" \
			InvokePipelineFunctionSchedule="$$(cat functions/environment.json | jq -r '.Functions.InvokePipeline.InvokePipelineFunctionSchedule')" \
			GitHubRepositoryOwner="$$(cat functions/environment.json | jq -r '.Functions.InvokePipeline.FunctionConfiguration.Environment.Variables.GITHUB_REPOSITORY_OWNER')" \
			GitHubRepositoryName="$$(cat functions/environment.json | jq -r '.Functions.InvokePipeline.FunctionConfiguration.Environment.Variables.GITHUB_REPOSITORY_NAME')" \
			PipelineStatePath="$$(cat functions/environment.json | jq -r '.Functions.InvokePipeline.FunctionConfiguration.Environment.Variables.PIPELINE_STATE_PATH')" \
			PipelineParamsPath="$$(cat functions/environment.json | jq -r '.Functions.InvokePipeline.FunctionConfiguration.Environment.Variables.PIPELINE_PARAMS_PATH')"  \
			InvokeLoadScriptFunctionName="$${STAGE}"-"$${APP_NAME}"-"$$(cat functions/environment.json | jq -r '.Functions.InvokeLoadScript.FunctionConfiguration.FunctionName')" \
			Stage="$${STAGE}" \
			AppName="$${APP_NAME}" \
			ECRBaseUri="${ECR_BASE_URI}" \
			BuildServiceRepositoryName="${BUILD_REPOSITORY}" \
			2>&1 | tee -a $$CFN_LOG_PATH || true

service.jobs.deploy:
	$(MAKE) -C jobs/ deploy

# TODO parameterize S3 config path and export as environment variable to recall in database shell scripts
service.config.deploy:
	@config_s3_path=s3://$$DATA_BUCKET_NAME/config/ && \
	echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying config to $$config_s3_path" 2>&1 | tee -a $$CFN_LOG_PATH && \
	aws s3 cp --recursive config/ $$config_s3_path 2>&1 | tee -a $$CFN_LOG_PATH

service.delete:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting ${SERVICE} service" 2>&1 | tee -a $$CFN_LOG_PATH
	$(MAKE) service.jobs.delete
	$(MAKE) service.functions.delete

service.functions.delete:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting CloudFormation" 2>&1 | tee -a $$CFN_LOG_PATH
	@aws cloudformation delete-stack \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" 2>&1 | tee -a $$CFN_LOG_PATH || true && \
	aws cloudformation wait stack-delete-complete \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" 2>&1 | tee -a $$CFN_LOG_PATH || true

service.jobs.delete:
	$(MAKE) -C jobs/ delete

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