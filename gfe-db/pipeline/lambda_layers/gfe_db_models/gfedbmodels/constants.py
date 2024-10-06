import os
import logging
from awsparameters import AppConfig
from awsparameters.manager import SessionManager

# Logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Environment variables
APP_NAME = os.environ["APP_NAME"]
STAGE = os.environ["STAGE"]
AWS_REGION = os.environ["AWS_REGION"]

session = SessionManager(region_name=AWS_REGION)
session.get_client("ssm")

# TODO parameterize paths or use a consolidated mapping for all layers
infra_config_path = f"/{APP_NAME}/{STAGE}/{AWS_REGION}/GfedbInfrastructureParamMappings"
infra = AppConfig(mappings_path=infra_config_path, boto3_session=session)

pipeline_config_path = f"/{APP_NAME}/{STAGE}/{AWS_REGION}/GfedbPipelineParamMappings"
pipeline = AppConfig(mappings_path=pipeline_config_path, boto3_session=session)

database_config_path = f"/{APP_NAME}/{STAGE}/{AWS_REGION}/GfedbDatabaseParamMappings"
database = AppConfig(mappings_path=database_config_path, boto3_session=session)
