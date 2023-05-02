"""
Loads gfe-db execution state to DynamoDB table.
"""
import os
import sys
sys.path.append(
    "/Users/ammon/Projects/nmdp-bioinformatics/02-Repositories/gfe-db/gfe-db/pipeline/functions/check_source_update"
)
from pathlib import Path
import logging
import boto3
from src.utils.types import (
    SourceConfig,
)
from src.utils.utils import (
    read_source_config,
    flatten_json_records,
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

ssm = boto3.client('ssm')
dynamodb = boto3.resource('dynamodb')

DATA_BUCKET_NAME = os.environ["DATA_BUCKET_NAME"]
PIPELINE_CONFIG_S3_PATH = os.environ["PIPELINE_CONFIG_S3_PATH"]
GITHUB_REPOSITORY_OWNER = os.environ["GITHUB_REPOSITORY_OWNER"]
GITHUB_REPOSITORY_NAME = os.environ["GITHUB_REPOSITORY_NAME"]

# fetch '/${AppName}/${Stage}/${AWS::Region}/ExecutionStateTableName' value from Parameter Store
execution_state_table_name = ssm.get_parameter(
    Name=f'/{os.environ["APP_NAME"]}/{os.environ["STAGE"]}/{os.environ["AWS_REGION"]}/ExecutionStateTableName'
)['Parameter']['Value']

# read in source config JSON file from local
source_config = read_source_config(DATA_BUCKET_NAME, PIPELINE_CONFIG_S3_PATH)
execution_history = source_config.repositories[f'{GITHUB_REPOSITORY_OWNER}/{GITHUB_REPOSITORY_NAME}'].execution_history

# flatten JSON records
execution_history_flat = flatten_json_records([item.dict() for item in execution_history])

# load to dynamodb table named execution_state_table_name using batch put
table = dynamodb.Table(execution_state_table_name)
with table.batch_writer() as batch:
    logger.info(f"Loading {len(execution_history_flat)} items to {execution_state_table_name}")
    for item in execution_history_flat:
        batch.put_item(Item=item)

logger.info(f"Success")