SERVICE := trigger

# # Local invoke
# invoke-post-job:
# 	sam local invoke \
# 		--event events/event.json \
# 		--env-vars local-env-vars.json PostJobGlueMetricsFunction \
# 		--profile ${PROFILE}

# Deploy
deploy:
	$(info [*] Packaging and deploying ${SERVICE} resources...)
	@sam build \
		--template-file template.yaml && \
	sam package \
		--resolve-s3 \
		--output-template-file packaged.yaml && \
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
			GitHubRepositoryName="$$(cat environment.json | jq -r '.FunctionConfiguration.Environment.Variables.GITHUB_REPOSITORY_NAME')"

delete: ##=> Delete service
	@aws cloudformation delete-stack \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" >/dev/null 2>&1
	@aws cloudformation wait stack-delete-complete \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" >/dev/null