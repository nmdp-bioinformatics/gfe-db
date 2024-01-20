"""In progress"""
import os

if __name__ != "app":
    import sys

    # for dev, local path to gfe-db modules
    # ./gfe-db/pipeline/lambda_layers/gfe_db_models (use absolute path)
    sys.path.append(os.environ["GFEDBMODELS_PATH"])
import logging
import json
import traceback
from gfedbmodels.types import ExecutionPayloadItem, ExecutionStateItem
from gfedbmodels.constants import session, pipeline
from gfedbmodels.utils import get_utc_now

logger = logging.getLogger()
logger.setLevel(logging.INFO)

dynamodb = session.resource("dynamodb")
table = dynamodb.Table(pipeline.params.GfeDbExecutionStateTableName)


def lambda_handler(event, context):
    logger.info(json.dumps(event))
    # return

    try:
        # validate input
        execution_payload_item = ExecutionPayloadItem(
            **json.loads(event["detail"]["input"])["input"]
        )
        status = event["detail"]["status"]

        # update execution state
        # composite key is commit.sha as commit__sha and execution.version as execution__version
        table.update_item(
            Key={
                "commit__sha": execution_payload_item.commit_sha,
                "execution__version": execution_payload_item.version,
            },
            UpdateExpression="SET #status = :status, #updated_utc = :updated_utc",
            ExpressionAttributeNames={
                "#status": "execution__status",
                "#updated_utc": "updated_utc",
            },
            ExpressionAttributeValues={
                ":status": status,
                ":updated_utc": get_utc_now(),
            },
        )

        return 0

    except Exception as e:
        logger.error(json.dumps(event))
        logger.error(traceback.format_exc())

        return 1


if __name__ == "__main__":
    from pathlib import Path

    event = json.loads((Path(__file__).parent / "event.json").read_text())
    lambda_handler(event, None)
