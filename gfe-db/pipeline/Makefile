SERVICE := pipeline
SCRIPTS_DIR := ${ROOT_DIR}/${APP_NAME}/${SERVICE}/scripts

target:
	$(info ${HELP_MESSAGE})
	@exit 0

# TODO: Don't deploy jobs if pipeline stack fails to create (exit Make)
service.deploy:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying ${SERVICE} service" 2>&1 | tee -a $$CFN_LOG_PATH
	$(MAKE) service.config.deploy
	$(MAKE) service.functions.deploy
	$(MAKE) service.jobs.deploy

service.functions.deploy:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying ${SERVICE} - functions" 2>&1 | tee -a $$CFN_LOG_PATH
	@sam build \
	--region "$${AWS_REGION}" \
		--use-container \
		--template-file template.yaml && \
	sam package \
	--region "$${AWS_REGION}" \
		--resolve-s3 \
		--output-template-file packaged.yaml && \
	sam deploy \
		--no-fail-on-empty-changeset \
		--region "$${AWS_REGION}" \
		--template-file packaged.yaml \
		--resolve-s3 \
		--stack-name "$${STAGE}-$${APP_NAME}-${SERVICE}" \
		--tags stage="$${STAGE}" app="$${APP_NAME}" service="${SERVICE}" branch="$$(git branch --show-current)" commit=$$(git rev-parse HEAD) \
		--capabilities CAPABILITY_IAM \
		--parameter-overrides \
			Stage="$${STAGE}" \
			AppName="$${APP_NAME}" \
			ServiceName="${SERVICE}" \
			createVpc="$${CREATE_VPC}" \
			usePrivateSubnet="$${USE_PRIVATE_SUBNET}" \
			ConfigS3Path="$${CONFIG_S3_PATH}" \
			GitHubRepositoryOwner="${GITHUB_REPOSITORY_OWNER}" \
			GitHubRepositoryName="${GITHUB_REPOSITORY_NAME}" \
			GitHubPersonalAccessToken="$$GITHUB_PERSONAL_ACCESS_TOKEN" \
			ECRBaseUri="${ECR_BASE_URI}" \
			BuildServiceRepositoryName="${BUILD_REPOSITORY_NAME}" \
			FeatureServiceUrl="${FEATURE_SERVICE_URL}" \
			Ec2KeyPairName="${EC2_KEY_PAIR_NAME}"

service.jobs.deploy:
	$(MAKE) -C jobs/ deploy

# TODO handle virtual environment creation and activation
service.state.build:
	@${PYTHON} ${SCRIPTS_DIR}/state/build.py ${ROOT_DIR}/${APP_NAME}/${SERVICE}/config/

service.state.load:
	@${PYTHON} ${SCRIPTS_DIR}/state/load.py ${ROOT_DIR}/${APP_NAME}/${SERVICE}/config/

# TODO parameterize S3 config path and export as environment variable to recall in database shell scripts
# TODO integrate and automate the build/load process for source config and execution state
service.config.deploy:
	$(MAKE) service.config.pipeline-params.deploy
	
service.config.pipeline-params.deploy:
	@config_s3_path=s3://${DATA_BUCKET_NAME}/${CONFIG_S3_PATH}/${SERVICE} && \
	echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying config to $$config_s3_path" 2>&1 | tee -a $$CFN_LOG_PATH && \
	aws s3 cp --recursive config/ $$config_s3_path 2>&1 | tee -a $$CFN_LOG_PATH

service.statemachine.update-pipeline.stop:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Stopping UpdatePipeline state machine" 2>&1 | tee -a $${CFN_LOG_PATH}
	@state_machine_arn=$$(aws ssm get-parameter \
		--name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/UpdatePipelineStateMachineArn" \
		--query "Parameter.Value" \
		--output text) && \
	executions=$$(aws stepfunctions list-executions \
		--state-machine-arn "$$state_machine_arn" \
		--status-filter RUNNING \
		--query "executions[*].executionArn" \
		--output text) && \
	if [ -n "$$executions" ]; then \
		for execution in $$executions; do \
			res=$$(aws stepfunctions stop-execution --execution-arn "$$execution") && \
			echo "$$res" | jq -r && \
			echo "\033[0;32mStopped execution:\033[0m" && \
			echo "\033[0;32m$$execution\033[0m"; \
		done; \
	else \
		echo "\033[0;32mNo running executions found for UpdatePipeline state machine\033[0m"; \
	fi
	
service.statemachine.load-concurrency-manager.stop:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Stopping LoadConcurrencyManager state machine" 2>&1 | tee -a $${CFN_LOG_PATH}
	@state_machine_arn=$$(aws ssm get-parameter \
		--name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/LoadConcurrencyManagerStateMachineArn" \
		--query "Parameter.Value" \
		--output text) && \
	executions=$$(aws stepfunctions list-executions \
		--state-machine-arn "$$state_machine_arn" \
		--status-filter RUNNING \
		--query "executions[*].executionArn" \
		--output text) && \
	if [ -n "$$executions" ]; then \
		for execution in $$executions; do \
			res=$$(aws stepfunctions stop-execution --execution-arn "$$execution") && \
			echo "$$res" | jq -r && \
			echo "\033[0;32mStopped execution:\033[0m" && \
			echo "\033[0;32m$$execution\033[0m"; \
		done; \
	else \
		echo "\033[0;32mNo running executions found for LoadConcurrencyManager state machine\033[0m"; \
	fi

service.queue.gfe-db-load.purge:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Purging GfeDbLoadQueue" 2>&1 | tee -a $${CFN_LOG_PATH}
	@queue_url=$$(aws ssm get-parameter \
		--name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/GfeDbLoadQueueUrl" \
		--query "Parameter.Value" \
		--output text) && \
	res=$$(aws sqs purge-queue --queue-url "$$queue_url") && \
	echo "$$res" | jq -r 2>&1 | tee -a $${CFN_LOG_PATH} && \
	echo "\033[0;32mGfeDbLoadQueue purged successfully\033[0m"

service.queue.gfe-db-processing.purge:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Purging GfeDbProcessingQueue" 2>&1 | tee -a $${CFN_LOG_PATH}
	@queue_url=$$(aws ssm get-parameter \
		--name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/GfeDbProcessingQueueUrl" \
		--query "Parameter.Value" \
		--output text) && \
	res=$$(aws sqs purge-queue --queue-url "$$queue_url") && \
	echo "$$res" | jq -r 2>&1 | tee -a $${CFN_LOG_PATH} && \
	echo "\033[0;32mGfeDbProcessingQueue purged successfully\033[0m"

service.alarm.update-pipeline-execution.status:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Fetching status of UpdatePipelineStateMachineExecutionAlarm" 2>&1 | tee -a $${CFN_LOG_PATH}
	@alarm_name=$$(aws ssm get-parameter \
		--name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/UpdatePipelineStateMachineExecutionAlarmName" \
		--query "Parameter.Value" \
		--output text) && \
	alarm_status=$$(aws cloudwatch describe-alarms \
		--alarm-names "$$alarm_name" \
		--query "MetricAlarms[0].StateValue" \
		--output text) && \
	echo "\033[0;32mUpdatePipelineStateMachineExecutionAlarm status: $$alarm_status\033[0m" 2>&1 | tee -a $${CFN_LOG_PATH}

service.alarm.update-pipeline-execution.wait: # arg: poll_interval (default: 10)
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Waiting for UpdatePipelineStateMachineExecutionAlarm status change" 2>&1 | tee -a $${CFN_LOG_PATH}
	@poll_interval=$${1:-10}; \
	alarm_name=$$(aws ssm get-parameter \
		--name "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/UpdatePipelineStateMachineExecutionAlarmName" \
		--query "Parameter.Value" \
		--output text) && \
	initial_status=$$(aws cloudwatch describe-alarms \
		--alarm-names "$$alarm_name" \
		--query "MetricAlarms[0].StateValue" \
		--output text) && \
	echo "Initial status: $$initial_status" && \
	while true; do \
		current_status=$$(aws cloudwatch describe-alarms \
			--alarm-names "$$alarm_name" \
			--query "MetricAlarms[0].StateValue" \
			--output text) && \
		if [ "$$current_status" != "$$initial_status" ]; then \
			echo "\033[0;32mAlarm status changed to: $$current_status\033[0m" 2>&1 | tee -a $${CFN_LOG_PATH} && \
			break; \
		fi; \
		echo "Current status: $$current_status. Waiting for $$poll_interval seconds..." && \
		sleep $$poll_interval; \
	done

service.delete:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting ${SERVICE} service" 2>&1 | tee -a $$CFN_LOG_PATH
	$(MAKE) service.jobs.delete
	$(MAKE) service.functions.delete

service.functions.delete:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting ${SERVICE} CloudFormation" 2>&1 | tee -a $$CFN_LOG_PATH
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