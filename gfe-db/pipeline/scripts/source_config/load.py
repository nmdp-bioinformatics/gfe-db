"""
Loads the initial gfe-db execution state to DynamoDB table.

TODO solution to avoid overwriting data when running this script (regular DynamoDB backups to S3 etc, fetch file from S3 and compare)
"""
import os
from pathlib import Path
import os
import sys
# for dev, local path to gfe-db modules
# ./gfe-db/pipeline/lambda_layers/gfe_db_models (use absolute path)
sys.path.append(os.environ["GFEDBMODELS_PATH"])

import logging
from datetime import datetime
utc_now = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
import json
import boto3
from gfedbmodels.constants import (
    execution_state_table_fields
)
from gfedbmodels.types import (
    ExecutionState,
)
from gfedbmodels.utils import (
    flatten_json_records,
    filter_null_fields
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

execution_state_table_name = ssm.get_parameter(
    Name=f'/{APP_NAME}/{STAGE}/{AWS_REGION}/ExecutionStateTableName'
)['Parameter']['Value']
table = dynamodb.Table(execution_state_table_name)

if __name__ == "__main__":

    # TODO scan table for existing items and throw error if not empty, require --overwrite flag to proceed

    # Paths
    input_dir = Path(sys.argv[1])

    # read in source config JSON file from local
    with open(input_dir / "execution-state.json", "r") as f:
        execution_state = ExecutionState(**json.load(f))

    # flatten JSON records for execution state table model
    # Uses double-underscore as separator because DynamoDB does not allow dots in attribute names
    execution_state_flat = flatten_json_records(
        execution_state.dict()["items"], 
        sep="__", 
        select_fields=[item.replace(".", "__") for item in execution_state_table_fields],
        filter_nulls=True)

    # load to dynamodb table named execution_state_table_name using batch put
    with table.batch_writer() as batch:
        logger.info(f"Loading {len(execution_state_flat)} items to {execution_state_table_name}")
        for item in execution_state_flat:
            batch.put_item(Item=item)

    logger.info(f"Success")