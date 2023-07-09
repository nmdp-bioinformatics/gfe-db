"""
Loads the initial gfe-db execution state to DynamoDB table.

TODO Sync state to local script so it can be reloaded from local
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
from gfedbmodels.constants import (
    session,
    pipeline
)
from gfedbmodels.types import (
    ExecutionState,
)
from gfedbmodels.utils import flatten_json_records

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

ssm = session.clients["ssm"]
dynamodb = session.resource("dynamodb")

# TODO
execution_state_table_fields = pipeline.params.ExecutionStateTableFields
execution_state_table_name = pipeline.params.ExecutionStateTableName

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
        select_fields=[
            item.replace(".", "__") for item in execution_state_table_fields
        ],
        filter_nulls=True,
    )

    # load to dynamodb table named execution_state_table_name using batch put
    with table.batch_writer() as batch:
        logger.info(
            f"Loading {len(execution_state_flat)} items to {execution_state_table_name}"
        )
        for item in execution_state_flat:
            batch.put_item(Item=item)

    logger.info(f"Success")
