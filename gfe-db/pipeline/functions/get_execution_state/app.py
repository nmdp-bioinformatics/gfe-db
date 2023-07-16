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

logger = logging.getLogger()
logger.setLevel(logging.INFO)

dynamodb = session.resource("dynamodb")
table = dynamodb.Table(pipeline.params.ExecutionStateTableName)

def lambda_handler(event, context):
    logger.info(json.dumps(event))

    # validate input
    execution_payload_item = ExecutionPayloadItem(**event)

    commit_state = table.get_item(
        Key={
            "commit__sha": execution_payload_item.commit_sha,
            "execution__version": execution_payload_item.version
        }
    )['Item']

    # Validate record with pydantic model
    execution_state_item = ExecutionStateItem.from_execution_state_item_json(commit_state)

    return execution_state_item.model_dump(exclude_none=True)


if __name__ == "__main__":
    from pathlib import Path

    event = json.loads((Path(__file__).parent / "event.json").read_text())

    lambda_handler(event, None)
