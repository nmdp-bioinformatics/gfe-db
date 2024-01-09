SERVICE := local

# Required environment variables
REQUIRED_VARS := STAGE APP_NAME AWS_REGION AWS_PROFILE DOCKER_USERNAME DOCKER_PASSWORD


env.validate.vars:
	$(foreach var,$(REQUIRED_VARS),\
		$(if $(value $(var)),,$(error $(var) is not set. Please add $(var) to the environment variables.)))

build: env.validate.vars
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Building ${SERVICE} ${APP_NAME}" 2>&1 | tee -a $$CFN_LOG_PATH
	$(MAKE) service.build

service.build:
	@from_path=$$(aws ssm get-parameter --name 	"/${APP_NAME}/${STAGE}/${AWS_REGION}/CurrentBackupS3Path" --query "Parameter.Value" --output text); \
	echo "Building ${SERVICE} from $$from_path" && \
	bash build-local.sh $$from_path
