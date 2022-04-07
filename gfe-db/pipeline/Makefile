SERVICE := pipeline

target:
	$(info ${HELP_MESSAGE})
	@exit 0

deploy:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying ${SERVICE} service" 2>&1 | tee -a $$CFN_LOG_PATH
	$(MAKE) deploy.pipeline 
	$(MAKE) deploy.jobs
	$(MAKE) deploy.config

deploy.pipeline:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying ${SERVICE} - functions" 2>&1 | tee -a $$CFN_LOG_PATH
	@sam build \
		--template-file template.yaml  2>&1 | tee -a $$CFN_LOG_PATH || true && \
	sam package \
		--resolve-s3 \
		--output-template-file packaged.yaml  2>&1 | tee -a $$CFN_LOG_PATH || true && \
	sam deploy \
		--no-fail-on-empty-changeset \
		--region "$${REGION}" \
		--template-file packaged.yaml \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" \
		--capabilities CAPABILITY_IAM \
		--parameter-overrides \
			Stage="$${STAGE}" \
			AppName="$${APP_NAME}" \
			TriggerFunctionName="$${STAGE}"-"$${APP_NAME}"-"$$(cat functions/environment.json | jq -r '.Functions.InvokePipeline.FunctionConfiguration.FunctionName')" \
			TriggerFunctionMemorySize="$$(cat functions/environment.json | jq -r '.Functions.InvokePipeline.FunctionConfiguration.MemorySize')" \
			TriggerFunctionTimeout="$$(cat functions/environment.json | jq -r '.Functions.InvokePipeline.FunctionConfiguration.Timeout')" \
			TriggerFunctionSchedule="$$(cat functions/environment.json | jq -r '.Functions.InvokePipeline.TriggerFunctionSchedule')" \
			GitHubRepositoryOwner="$$(cat functions/environment.json | jq -r '.Functions.InvokePipeline.FunctionConfiguration.Environment.Variables.GITHUB_REPOSITORY_OWNER')" \
			GitHubRepositoryName="$$(cat functions/environment.json | jq -r '.Functions.InvokePipeline.FunctionConfiguration.Environment.Variables.GITHUB_REPOSITORY_NAME')" \
			PipelineStatePath="$$(cat functions/environment.json | jq -r '.Functions.InvokePipeline.FunctionConfiguration.Environment.Variables.PIPELINE_STATE_PATH')" \
			PipelineParamsPath="$$(cat functions/environment.json | jq -r '.Functions.InvokePipeline.FunctionConfiguration.Environment.Variables.PIPELINE_PARAMS_PATH')"  \
			LoadFunctionName="$${STAGE}"-"$${APP_NAME}"-"$$(cat functions/environment.json | jq -r '.Functions.InvokeLoadScript.FunctionConfiguration.FunctionName')" \
			Stage="$${STAGE}" \
			AppName="$${APP_NAME}" \
			ECRBaseUri="${ECR_BASE_URI}" \
			BuildServiceRepositoryName="${BUILD_REPOSITORY}" \
			LoadServiceRepositoryName="${LOAD_REPOSITORY}" \
			2>&1 | tee -a $$CFN_LOG_PATH || true

# deploy.cfn:
# 	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying CloudFormation" 2>&1 | tee -a $$CFN_LOG_PATH
# 	@aws cloudformation deploy \
# 		--region $${REGION} \
# 		--template-file template.yaml \
# 		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" \
# 		--capabilities CAPABILITY_IAM \
# 		--parameter-overrides \
# 			Stage="$${STAGE}" \
# 			AppName="$${APP_NAME}" \
# 			ECRBaseUri="${ECR_BASE_URI}" \
# 			BuildServiceRepositoryName="${BUILD_REPOSITORY}" \
# 			LoadServiceRepositoryName="${LOAD_REPOSITORY}" 2>&1 | tee -a $$CFN_LOG_PATH

deploy.jobs:
	$(MAKE) -C jobs/ deploy

# deploy.functions:
# 	$(MAKE) -C functions/ deploy

deploy.config:
	@config_s3_path=s3://$$DATA_BUCKET_NAME/config/ && \
	echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying config to $$config_s3_path" 2>&1 | tee -a $$CFN_LOG_PATH && \
	aws s3 cp --recursive config/ $$config_s3_path 2>&1 | tee -a $$CFN_LOG_PATH

delete:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting ${SERVICE} service" 2>&1 | tee -a $$CFN_LOG_PATH
	$(MAKE) delete.jobs
	$(MAKE) delete.pipeline

delete.pipeline:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting CloudFormation" 2>&1 | tee -a $$CFN_LOG_PATH
	@aws cloudformation delete-stack \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" 2>&1 | tee -a $$CFN_LOG_PATH || true && \
	aws cloudformation wait stack-delete-complete \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" 2>&1 | tee -a $$CFN_LOG_PATH || true

delete.jobs:
	$(MAKE) -C jobs/ delete

# delete.functions:
# 	$(MAKE) -C functions/ delete

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