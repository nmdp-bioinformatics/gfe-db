SERVICE := pipeline

target:
	$(info ${HELP_MESSAGE})
	@exit 0

# TODO: Don't deploy jobs if pipeline stack fails to create (exit Make)
deploy:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying ${SERVICE} service" 2>&1 | tee -a $$CFN_LOG_PATH
	$(MAKE) service.config.deploy
	$(MAKE) service.deploy
	$(MAKE) service.jobs.deploy

service.deploy:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying ${SERVICE} - functions" 2>&1 | tee -a $$CFN_LOG_PATH
	@sam build \
		--use-container \
		--template-file template.yaml && \
	sam package \
		--resolve-s3 \
		--output-template-file packaged.yaml && \
	sam deploy \
		--no-fail-on-empty-changeset \
		--region "$${AWS_REGION}" \
		--template-file packaged.yaml \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" \
		--tags stage="$${STAGE}" app="$${APP_NAME}" service="${SERVICE}" branch="$$(git branch --show-current)" commit=$$(git rev-parse HEAD) \
		--capabilities CAPABILITY_IAM \
		--parameter-overrides \
			Stage="$${STAGE}" \
			AppName="$${APP_NAME}" \
			createVpc="$${CREATE_VPC}" \
			usePrivateSubnet="$${USE_PRIVATE_SUBNET}" \
			GitHubRepositoryOwner="${GITHUB_REPOSITORY_OWNER}" \
			GitHubRepositoryName="${GITHUB_REPOSITORY_NAME}" \
			InvokePipelineFunctionName="$${STAGE}"-"$${APP_NAME}"-"$$(cat functions/environment.json | jq -r '.Functions.InvokePipeline.FunctionConfiguration.FunctionName')" \
			InvokePipelineFunctionSchedule="$$(cat functions/environment.json | jq -r '.Functions.InvokePipeline.InvokePipelineFunctionSchedule')" \
			InvokePipelineFunctionMemorySize="$$(cat functions/environment.json | jq -r '.Functions.InvokePipeline.FunctionConfiguration.MemorySize')" \
			InvokePipelineFunctionTimeout="$$(cat functions/environment.json | jq -r '.Functions.InvokePipeline.FunctionConfiguration.Timeout')" \
			PipelineStatePath="$$(cat functions/environment.json | jq -r '.Functions.InvokePipeline.FunctionConfiguration.Environment.Variables.PIPELINE_STATE_PATH')" \
			PipelineParamsPath="$$(cat functions/environment.json | jq -r '.Functions.InvokePipeline.FunctionConfiguration.Environment.Variables.PIPELINE_PARAMS_PATH')" \
			ExecuteValidationQueriesFunctionName="$${STAGE}"-"$${APP_NAME}"-"$$(cat functions/environment.json | jq -r '.Functions.ExecuteValidationQueries.FunctionConfiguration.FunctionName')" \
			ExecuteValidationQueriesFunctionMemorySize="$$(cat functions/environment.json | jq -r '.Functions.ExecuteValidationQueries.FunctionConfiguration.MemorySize')" \
			ExecuteValidationQueriesFunctionTimeout="$$(cat functions/environment.json | jq -r '.Functions.ExecuteValidationQueries.FunctionConfiguration.Timeout')" \
			InvokeBackupScriptFunctionName="$${STAGE}"-"$${APP_NAME}"-"$$(cat functions/environment.json | jq -r '.Functions.InvokeBackupScript.FunctionConfiguration.FunctionName')" \
			InvokeLoadScriptFunctionName="$${STAGE}"-"$${APP_NAME}"-"$$(cat functions/environment.json | jq -r '.Functions.InvokeLoadScript.FunctionConfiguration.FunctionName')" \
			ValidateBuildOutputFunctionName="$${STAGE}"-"$${APP_NAME}"-"$$(cat functions/environment.json | jq -r '.Functions.ValidateBuildOutput.FunctionConfiguration.FunctionName')" \
			DisableBackupFunctionName="$${STAGE}"-"$${APP_NAME}"-"$$(cat functions/environment.json | jq -r '.Functions.DisableBackup.FunctionConfiguration.FunctionName')" \
			ECRBaseUri="${ECR_BASE_URI}" \
			BuildServiceRepositoryName="${BUILD_REPOSITORY_NAME}" \
			FeatureServiceUrl="${FEATURE_SERVICE_URL}" \
			2>&1 | tee -a $$CFN_LOG_PATH || true

service.jobs.deploy:
	$(MAKE) -C jobs/ deploy

service.config.deploy:
	$(MAKE) service.config.pipeline-params.deploy

service.config.pipeline-params.deploy:
	@config_s3_path=s3://$$DATA_BUCKET_NAME/config/pipeline/ && \
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