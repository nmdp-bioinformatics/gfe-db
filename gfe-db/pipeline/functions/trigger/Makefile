SERVICE := pipeline-trigger

# # Local invoke
# invoke-post-job:
# 	sam local invoke \
# 		--event events/event.json \
# 		--env-vars local-env-vars.json PostJobGlueMetricsFunction \
# 		--profile ${PROFILE}

# Deploy
deploy:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying ${SERVICE}" 2>&1 | tee -a $$CFN_LOG_PATH
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
			TriggerFunctionName="$${STAGE}"-"$${APP_NAME}"-"${SERVICE}" \
			TriggerFunctionMemorySize="$$(cat environment.json | jq -r '.FunctionConfiguration.MemorySize')" \
			TriggerFunctionTimeout="$$(cat environment.json | jq -r '.FunctionConfiguration.Timeout')" \
			TriggerFunctionSchedule="$$(cat environment.json | jq -r '.TriggerFunctionSchedule')" \
			GitHubRepositoryOwner="$$(cat environment.json | jq -r '.FunctionConfiguration.Environment.Variables.GITHUB_REPOSITORY_OWNER')" \
			GitHubRepositoryName="$$(cat environment.json | jq -r '.FunctionConfiguration.Environment.Variables.GITHUB_REPOSITORY_NAME')" \
			PipelineStatePath="$$(cat environment.json | jq -r '.FunctionConfiguration.Environment.Variables.PIPELINE_STATE_PATH')" \
			PipelineParamsPath="$$(cat environment.json | jq -r '.FunctionConfiguration.Environment.Variables.PIPELINE_PARAMS_PATH')"  \
			2>&1 | tee -a $$CFN_LOG_PATH || true

delete: ##=> Delete service
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting ${SERVICE}" 2>&1 | tee -a $$CFN_LOG_PATH
	@aws cloudformation delete-stack \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" 2>&1 | tee -a $$CFN_LOG_PATH || true
	@aws cloudformation wait stack-delete-complete \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" 2>&1 | tee -a $$CFN_LOG_PATH || true