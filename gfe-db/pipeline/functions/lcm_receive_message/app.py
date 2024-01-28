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
    database)

# set up logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Environment
APP_NAME = os.environ["APP_NAME"]
STAGE = os.environ["STAGE"]
AWS_REGION = os.environ["AWS_REGION"]

# Boto3 Clients
sqs = session.client("sqs", region_name=AWS_REGION)

# Get SSM Parameters
gfe_db_load_queue_url = database.params.GfeDbLoadQueueUrl

def lambda_handler(event, context):

    logger.info(json.dumps(event))

    response = sqs.receive_message(
        QueueUrl=gfe_db_load_queue_url,
        MaxNumberOfMessages=1
    )

    if "Messages" in response:
        message = response["Messages"][0]
    else:
        return_msg = "No messages found in GfeDbLoadQueue."
        logger.info(return_msg)
        return {}
    
    return message


if __name__ == "__main__":
    from pathlib import Path

    event = json.loads((Path(__file__).parent / "sqs-event.json").read_text())
    lambda_handler(event, "")
