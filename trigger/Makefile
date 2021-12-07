SERVICE := trigger
STAGE := dev
APP_NAME := gfe-db
REGION := us-east-1

# # Local invoke
# invoke-post-job:
# 	sam local invoke \
# 		--event events/event.json \
# 		--env-vars local-env-vars.json PostJobGlueMetricsFunction \
# 		--profile ${PROFILE}

# Deploy
deploy: ##==> command line args: environment=dev,qa,prod etc; example: make deploy environment=qa
	$(info [*] Packaging and deploying ${SERVICE} resources...)
	@rm -rf .aws-sam
	@sam build && \
	sam deploy \
		--no-fail-on-empty-changeset \
		--resolve-s3 \
		--region "${REGION}" \
		--template-file template.yaml \
		--stack-name "${APP_NAME}-${SERVICE}" \
		--capabilities CAPABILITY_IAM \
		--parameter-overrides \
			Stage="${STAGE}" \
			AppName="${APP_NAME}" \
			TriggerFunctionName="${STAGE}"-"${APP_NAME}"-"${SERVICE}" \
			GitHubPersonalAccessToken="$$(cat environment.json | jq -r '.FunctionConfiguration.Environment.Variables.GITHUB_PERSONAL_ACCESS_TOKEN')" \
			GitHubRepositoryOwner="$$(cat environment.json | jq -r '.FunctionConfiguration.Environment.Variables.GITHUB_REPOSITORY_OWNER')" \
			GitHubRepositoryName="$$(cat environment.json | jq -r '.FunctionConfiguration.Environment.Variables.GITHUB_REPOSITORY_NAME')" \
			DataBucketName="$$(cat environment.json | jq -r '.FunctionConfiguration.Environment.Variables.GFE_BUCKET')"

delete: ##=> Delete service
	sam delete --stack-name ${APP_NAME}-${SERVICE} >/dev/null