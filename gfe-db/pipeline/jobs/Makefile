target:
	$(info ${HELP_MESSAGE})
	@exit 0

deploy: ecr.login deploy.jobs.build deploy.jobs.load

ecr.login:
	$(info [*] Logging into ECR...)
	@aws ecr get-login-password \
		--region "$${REGION}" | docker login \
			--username AWS \
			--password-stdin "$${ECR_BASE_URI}"

deploy.jobs.build:
	$(info [*] Pushing Build Service Docker image to ECR...)
	@docker build -t "$${BUILD_REPOSITORY}" build/ && \
	docker tag "$${BUILD_REPOSITORY}:latest" "$${ECR_BASE_URI}/$${BUILD_REPOSITORY}:latest" && \
	docker push "$${ECR_BASE_URI}/$${BUILD_REPOSITORY}:latest"

deploy.jobs.load:
	$(info [*] Pushing Build Service Docker image to ECR...)
	@docker build -t "$${LOAD_REPOSITORY}" load/ && \
	docker tag "$${LOAD_REPOSITORY}:latest" "$${ECR_BASE_URI}/$${LOAD_REPOSITORY}:latest" && \
	docker push "$${ECR_BASE_URI}/$${LOAD_REPOSITORY}:latest"

delete: delete.jobs.build delete.jobs.load

delete.jobs.build:
	$(info [*] Deleting Build Service Docker images...)
	@ecr_images="$$(aws ecr list-images --region $${REGION} --repository-name $${BUILD_REPOSITORY} --query 'imageIds[*]' --output json)" && \
	[ "$$ecr_images" = "[]" ] && \
	echo "No images to delete in $${ECR_BASE_URI}/$${BUILD_REPOSITORY}" && exit 0 || \
	echo "Deleting images in $${ECR_BASE_URI}/$${BUILD_REPOSITORY}" && \
	aws ecr batch-delete-image \
		--region "$${REGION}" \
		--repository-name "$${BUILD_REPOSITORY}" \
		--image-ids "$$ecr_images" || true

delete.jobs.load:
	$(info [*] Deleting Load Service Docker images...)
	@ecr_images="$$(aws ecr list-images --region $${REGION} --repository-name $${LOAD_REPOSITORY} --query 'imageIds[*]' --output json)" && \
	[ "$$ecr_images" = "[]" ] && \
	echo "No images to delete in $${ECR_BASE_URI}/$${LOAD_REPOSITORY}" && exit 0 || \
	echo "Deleting images in $${ECR_BASE_URI}/$${LOAD_REPOSITORY}" && \
	aws ecr batch-delete-image \
		--region "$${REGION}" \
		--repository-name "$${LOAD_REPOSITORY}" \
		--image-ids "$$ecr_images" || true

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