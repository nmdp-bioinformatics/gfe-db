import os
import logging
import boto3
from awsparameters import AppConfig

# Logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Environment variables
APP_NAME = os.environ["APP_NAME"]
STAGE = os.environ["STAGE"]
AWS_REGION = os.environ["AWS_REGION"]

session = boto3.Session(region_name=AWS_REGION)

infra_config_path = f"/{APP_NAME}/{STAGE}/{AWS_REGION}/GfedbInfrastructureParamMappings"
infra = AppConfig(mapping_path=infra_config_path, boto3_session=session)

pipeline_config_path = f"/{APP_NAME}/{STAGE}/{AWS_REGION}/GfedbPipelineParamMappings"
pipeline = AppConfig(mapping_path=pipeline_config_path, boto3_session=session)

GITHUB_PERSONAL_ACCESS_TOKEN = pipeline.secrets.git_hub_personal_access_token
github_source_repository = pipeline.params.git_hub_source_repository
GITHUB_REPOSITORY_OWNER = pipeline.params.git_hub_source_repository["owner"]
GITHUB_REPOSITORY_NAME = pipeline.params.git_hub_source_repository["name"]
execution_state_table_name = pipeline.params.execution_state_table_name
# state_machine_arn = pipeline.params.state_machine_arn
data_bucket_name = infra.params.data_bucket_name
gfedb_processing_queue_url = pipeline.params.gfe_db_processing_queue_url
execution_state_table_fields = pipeline.params.execution_state_table_fields
