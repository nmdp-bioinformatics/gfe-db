# PROFILE := 531868584498_AdministratorAccess
# SERVICE := build

# # Local invoke
# invoke-getConfigTableData:
# 	sam local invoke \
# 		--event src/getConfigTableData/event.json \
# 		--env-vars local-env-vars.json GetConfigTableDataFunction \
# 		--profile ${PROFILE}

# # Deploy
# deploy: ##=> Deploy resources using SAM
# 	$(info [*] Packaging and deploying ${SERVICE} resources...)
# 	# sam build && \ # Needed for Lambda layers
# 	sam package \
# 		--s3-bucket $${DEPLOYMENT_BUCKET_NAME} \
# 		--output-template-file packaged.yaml && \
# 	sam deploy \
# 		--no-fail-on-empty-changeset \
# 		--region $${REGION} \
# 		--template-file packaged.yaml \
# 		--stack-name $${STACK_NAME}-${SERVICE}-$${AWS_BRANCH} \
# 		--capabilities CAPABILITY_IAM \
# 		--parameter-overrides \
# 			Stage=$${AWS_BRANCH} \
# 			AppName=$${STACK_NAME} \
# 			RestApiName=$$(cat local-env-vars.json | jq -r '.RestApi.RestApiName')

# delete: ##=> Delete service
# 	aws cloudformation delete-stack --stack-name $${STACK_NAME}-${SERVICE}-$${AWS_BRANCH} >/dev/null