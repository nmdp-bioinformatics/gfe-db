"""
Loads the initial gfe-db execution state to DynamoDB table.
"""
import os
import sys
sys.path.append(
    "/Users/ammon/Projects/nmdp-bioinformatics/02-Repositories/gfe-db/gfe-db/pipeline/functions/check_source_update"
)
from pathlib import Path
import logging
from datetime import datetime
utc_now = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
import json
import boto3
from src.utils.constants import (
    execution_state_table_fields
)
from src.utils.types import (
    ExecutionState,
)
from src.utils.utils import (
    flatten_json_records
)

logging.basicConfig(level=logging.INFO)
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

    # TODO scan table for existing items and throw error if not empty

    # Paths
    # output_dir = Path(f"{APP_NAME}/pipeline/config")
    output_dir = Path(sys.argv[1])

    # read in source config JSON file from local
    with open(output_dir / "execution-state.json", "r") as f:
        execution_state = ExecutionState(**json.load(f))

    # execution_state_json = [
    #     {
    #         "commit.sha": item.commit.sha,
    #         "version": item.execution.version,
    #         "execution": item.execution.json(),
    #         "commit": item.commit.json(),
    #         "repository": item.repository.json(),
    #     } for item in execution_state.items
    # ]

    # TODO use selected fields from constants
    # flatten JSON records and filter nulls
    # skip_fields = [
    #     "execution.input_parameters",
    #     "repository.description",
    #     "repository.excluded_commit_shas",
    #     "repository.target_metadata_config",
    #     "repository.tracked_assets"

    # ]
    execution_state_flat = flatten_json_records(execution_state.dict()["items"], sep="__", select_fields=[item.replace(".", "__") for item in execution_state_table_fields])

    # load to dynamodb table named execution_state_table_name using batch put
    table = dynamodb.Table(execution_state_table_name)
    with table.batch_writer() as batch:
        logger.info(f"Loading {len(execution_state_flat)} items to {execution_state_table_name}")
        for item in execution_state_flat:
            batch.put_item(Item=item)

    logger.info(f"Success")