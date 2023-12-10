SERVICE := pipeline-jobs

target:
	$(info ${HELP_MESSAGE})
	@exit 0

deploy:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying ${SERVICE}" 2>&1 | tee -a $$CFN_LOG_PATH
	$(MAKE) deploy.ecr.jobs.build

deploy.ecr.jobs.build:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Pushing Build Service Docker image to ECR" 2>&1 | tee -a $$CFN_LOG_PATH
	@aws ecr get-login-password \
		--profile "$${AWS_PROFILE}" \
		--region "$${AWS_REGION}" | docker login \
			--username AWS \
			--password-stdin "$${ECR_BASE_URI}" 2>&1 | tee -a $$CFN_LOG_PATH || true && \
	docker build \
		-t "$${BUILD_REPOSITORY_NAME}" \
		--platform "linux/amd64" \
		build/ && \
	docker tag "$${BUILD_REPOSITORY_NAME}:latest" "$${ECR_BASE_URI}/$${BUILD_REPOSITORY_NAME}:latest" && \
	docker push "$${ECR_BASE_URI}/$${BUILD_REPOSITORY_NAME}:latest" 2>&1 | tee -a $$CFN_LOG_PATH || true

delete:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting ${SERVICE}" 2>&1 | tee -a $$CFN_LOG_PATH
	$(MAKE) delete.ecr.jobs.build

delete.ecr.jobs.build:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting Build Service Docker images" 2>&1 | tee -a $$CFN_LOG_PATH
	@ecr_images="$$(aws ecr list-images \
		--region $${AWS_REGION} \
		--repository-name $${BUILD_REPOSITORY_NAME} \
		--query 'imageIds[*]' \
		--output json || true)" && \
	[ "$$ecr_images" = "[]" ] && \
	echo "No images to delete in $${ECR_BASE_URI}/$${BUILD_REPOSITORY_NAME}" && exit 0 || \
	echo "Deleting images in $${ECR_BASE_URI}/$${BUILD_REPOSITORY_NAME}" && \
	aws ecr batch-delete-image \
		--region "$${AWS_REGION}" \
		--repository-name "$${BUILD_REPOSITORY_NAME}" \
		--image-ids "$$ecr_images" 2>&1 | tee -a $$CFN_LOG_PATH || true

#############
#  Helpers  #
#############

define HELP_MESSAGE

	Common usage:

	...::: Deploy all CloudFormation based services :::...
	$ make deploy

	...::: Delete all CloudFormation based services :::...
	$ make delete

endef