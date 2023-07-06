import os
import boto3
from awsparameters import AppConfig
from gfedbmodels.constants import (
    session,
    infra,
    pipeline
)

# Environment
AWS_REGION = os.environ["AWS_REGION"]
APP_NAME = os.environ["APP_NAME"]
STAGE = os.environ["STAGE"]
PIPELINE_SOURCE_CONFIG_S3_PATH = os.environ["PIPELINE_SOURCE_CONFIG_S3_PATH"]

(
    GITHUB_REPOSITORY_OWNER, 
    GITHUB_REPOSITORY_NAME, 
    execution_state_table_name, 
    gfedb_processing_queue_url, 
    execution_state_table_fields
) = (
    pipeline.params.GitHubSourceRepository["owner"],
    pipeline.params.GitHubSourceRepository["name"],
    pipeline.params.ExecutionStateTableName,
    pipeline.params.GfeDbProcessingQueueUrl,
    pipeline.params.ExecutionStateTableFields
)

data_bucket_name = infra.params.DataBucketName