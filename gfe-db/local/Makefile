SERVICE := local

build:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Building ${SERVICE} ${APP_NAME}" 2>&1 | tee -a $$CFN_LOG_PATH
	$(MAKE) service.build

service.build:
	@from_path=$$(aws ssm get-parameter --name 	"/${APP_NAME}/${STAGE}/${AWS_REGION}/CurrentBackupS3Path" --query "Parameter.Value" --output text); \
	echo "Building ${SERVICE} from $$from_path" && \
	bash build-local.sh $$from_path
