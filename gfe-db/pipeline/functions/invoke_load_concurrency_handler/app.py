"""This function is invoked through SNS when GfeDbLoadQueueHasMessagesAlarm is triggered from messages present in the GfeDbLoadQueue.
It's only responsibility is to invoke the LoadConcurrencyHandler state machine which maintains a concurrency of 1 for loading Neo4j
to avoid clashes with concurrent executions in GfeDbUpdatePipeline. This allows GfeDbUpdatePipeline to run data builds at concurrency > 1
and keeps data loads at concurrency = 1. The LoadConcurrencyHandler will end the execution when GfeDbLoadQueueHasMessagesAlarm enters 
OK state, meaning there are no more messages in the queue. All requests for loading Neo4j are handled by the LoadConcurrencyHandler 
state machine.
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

# set up logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Environment
APP_NAME = os.environ["APP_NAME"]
STAGE = os.environ["STAGE"]
AWS_REGION = os.environ["AWS_REGION"]

# Boto3 Clients
states = session.client("stepfunctions", region_name=AWS_REGION)

# # TODO Get SSM Parameters
# lch_state_machine_arn = pipeline.params.LoadConcurrencyHandlerStateMachineArn

def lambda_handler(event, context):

    logger.info(json.dumps(event))

    alarm_message = json.loads(event["Records"][0]["Sns"]["Message"])

    # validate the alarm state is IN ALARM
    state_has_changed = "NewStateValue" in alarm_message
    is_in_alarm = alarm_message["NewStateValue"] == "ALARM"

    load_queue_has_messages = state_has_changed and is_in_alarm

    if load_queue_has_messages:

        # TODO query the state table for commits with PENDING status

        # TODO trigger the Load Concurrency Handler state machine
        pass

    return
    # return {
    #     "statusCode": 200,
    #     "body": json.dumps({"message": return_msg, "execution_arns": execution_arns}),
    # }


if __name__ == "__main__":
    from pathlib import Path

    event = json.loads((Path(__file__).parent / "sns-event.json").read_text())
    lambda_handler(event, "")
