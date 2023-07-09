"""In progress"""
import os

if __name__ != "app":
    import sys

    # for dev, local path to gfe-db modules
    # ./gfe-db/pipeline/lambda_layers/gfe_db_models (use absolute path)
    sys.path.append(os.environ["GFEDBMODELS_PATH"])
import logging
import json
from gfedbmodels.types import (
    ExecutionPayloadItem,
    ExecutionStateItem
)
from gfedbmodels.constants import (
    session,
    pipeline
)
from gfedbmodels.utils import (
    restore_nested_json
)

logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Environment
APP_NAME = os.environ["APP_NAME"]
STAGE = os.environ["STAGE"]

dynamodb = session.resource("dynamodb")
table = dynamodb.Table(pipeline.params.ExecutionStateTableName)

def lambda_handler(event, context):
    logger.info(json.dumps(event))

    # validate input
    execution_payload_item = ExecutionPayloadItem(**event)

    # TODO Get state for commit in input
    commit_state = table.get_item(
        Key={
            "commit__sha": execution_payload_item.commit_sha,
            "execution__version": execution_payload_item.version
        }
    )['Item']

    commit_state = restore_nested_json(commit_state, split_on="__")
    commit_state = ExecutionStateItem(**commit_state) # TODO table logic in models.utils

    # TODO Return state, include the SQS message receipt in case it needs to be returned to the queue if the state machine fails

    return




if __name__ == "__main__":
    from pathlib import Path

    event = json.loads((Path(__file__).parent / "event.json").read_text())

    lambda_handler(event, None)
