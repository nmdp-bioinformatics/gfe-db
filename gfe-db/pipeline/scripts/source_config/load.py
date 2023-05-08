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
import json
import boto3
from src.utils.types import (
    ExecutionState,
)
from src.utils.utils import (
    flatten_json_records
)

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

ssm = boto3.client('ssm')
dynamodb = boto3.resource('dynamodb')

AWS_REGION = os.environ["AWS_REGION"]
APP_NAME = os.environ["APP_NAME"]
STAGE = os.environ["STAGE"]
DATA_BUCKET_NAME = os.environ["DATA_BUCKET_NAME"]
GITHUB_REPOSITORY_OWNER = os.environ["GITHUB_REPOSITORY_OWNER"]
GITHUB_REPOSITORY_NAME = os.environ["GITHUB_REPOSITORY_NAME"]

# fetch '/${AppName}/${Stage}/${AWS::Region}/ExecutionStateTableName' value from Parameter Store
execution_state_table_name = ssm.get_parameter(
    Name=f'/{APP_NAME}/{STAGE}/{AWS_REGION}/ExecutionStateTableName'
)['Parameter']['Value']

if __name__ == "__main__":

    # Paths
    # TODO arg
    output_dir = Path(f"{APP_NAME}/pipeline/config")

    # read in source config JSON file from local
    with open(output_dir / "execution-state.json", "r") as f:
        execution_state = ExecutionState(**json.load(f))

    # TODO BOOKMARK redeploy execution state table
    # flatten JSON records and filter nulls
    execution_state_flat = [
        {
            k: v for k, v in item.items() if v is not None
        } for item in flatten_json_records([item.dict() for item in execution_state.items])
    ]

    # load to dynamodb table named execution_state_table_name using batch put
    table = dynamodb.Table(execution_state_table_name)
    with table.batch_writer() as batch:
        logger.info(f"Loading {len(execution_state_flat)} items to {execution_state_table_name}")
        for item in execution_state_flat:
            batch.put_item(Item=item)

    logger.info(f"Success")