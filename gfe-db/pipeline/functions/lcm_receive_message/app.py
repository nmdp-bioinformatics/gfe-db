"""This function is invoked through by the LoadConcurrencyManager state machine when it is triggered by the GfeDbLoadQueueHasMessagesAlarm.
It polls the GfeDbLoadQueue for messages and invokes the LoadNeo4j state machine for each message. If no messages are found the state machine will check the alarm status
and repeat the polling process until the alarm is in OK state.
"""

import os
if __name__ != "app":
    import sys

    # for dev, local path to gfe-db modules
    # ./gfe-db/pipeline/lambda_layers/gfe_db_models (use absolute path)
    sys.path.append(os.environ["GFEDBMODELS_PATH"])

import logging
import json
from gfedbmodels.constants import (
    session,
    pipeline)

logger = logging.getLogger()
logger.setLevel(logging.INFO)

APP_NAME = os.environ["APP_NAME"]
STAGE = os.environ["STAGE"]
AWS_REGION = os.environ["AWS_REGION"]

sqs = session.client("sqs", region_name=AWS_REGION)

gfe_db_load_queue_url = pipeline.params.GfeDbLoadQueueUrl

def lambda_handler(event, context):

    logger.info(json.dumps(event))

    res = sqs.receive_message(
        QueueUrl=gfe_db_load_queue_url,
        MaxNumberOfMessages=1
    )

    if "Messages" in res:
        message = res["Messages"][0]

        # Format the message body as json
        message['Body'] = json.loads(message['Body'])

        # change message visibility to 8 hours
        sqs.change_message_visibility(
            QueueUrl=gfe_db_load_queue_url,
            ReceiptHandle=message["ReceiptHandle"],
            VisibilityTimeout=28800
        )

    else:
        logger.info("No messages found")
        return {}
    
    return message


if __name__ == "__main__":
    from pathlib import Path

    event = json.loads((Path(__file__).parent / "event.json").read_text())
    lambda_handler(event, "")
