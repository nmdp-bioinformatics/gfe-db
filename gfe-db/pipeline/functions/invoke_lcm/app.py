"""This function is invoked through SNS when GfeDbLoadQueueHasMessagesAlarm is triggered from messages present in the GfeDbLoadQueue.
It's only responsibility is to invoke the LoadConcurrencyManager state machine which maintains a concurrency of 1 for loading Neo4j
to avoid clashes with concurrent executions in GfeDbUpdatePipeline. This allows GfeDbUpdatePipeline to run data builds at concurrency > 1
and keeps data loads at concurrency = 1. The LoadConcurrencyManager will end the execution when GfeDbLoadQueueHasMessagesAlarm enters 
OK state, meaning there are no more messages in the queue. All requests for loading Neo4j are handled by the LoadConcurrencyManager 
state machine.
"""

import os
if __name__ != "app":
    import sys

    # for dev, local path to gfe-db modules
    # ./gfe-db/pipeline/lambda_layers/gfe_db_models (use absolute path)
    sys.path.append(os.environ["GFEDBMODELS_PATH"])

import logging
from datetime import datetime
import json
from gfedbmodels.constants import (
    session,
    pipeline)

logger = logging.getLogger()
logger.setLevel(logging.INFO)

APP_NAME = os.environ["APP_NAME"]
STAGE = os.environ["STAGE"]
AWS_REGION = os.environ["AWS_REGION"]

states = session.client("stepfunctions", region_name=AWS_REGION)

update_pipeline_state_machine_arn = pipeline.params.UpdatePipelineStateMachineArn
lcm_state_machine_arn = pipeline.params.LoadConcurrencyManagerStateMachineArn

def lambda_handler(event, context):

    logger.info(json.dumps(event))

    alarm_message = json.loads(event["Records"][0]["Sns"]["Message"])

    # Validate the alarm state is IN ALARM
    state_has_changed = "NewStateValue" in alarm_message
    is_in_alarm = alarm_message["NewStateValue"] == "ALARM"

    if state_has_changed and is_in_alarm:

        # TODO query the state table for commits with PENDING status to get the invocation_id for the LCM's execution_id
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        execution_id = update_pipeline_state_machine_arn.split(":")[-1] + "_" + timestamp

        if not executions_in_progress(lcm_state_machine_arn):
            response = states.start_execution(
                stateMachineArn=lcm_state_machine_arn,
                name=execution_id
            )

    return {
        "statusCode": 200
    }


def executions_in_progress(state_machine_arn):
    # List executions for the state machine
    response = states.list_executions(
        stateMachineArn=state_machine_arn,
        statusFilter="RUNNING"
    )

    return bool(response['executions'])


if __name__ == "__main__":
    from pathlib import Path

    event = json.loads((Path(__file__).parent / "gfedbloadqueue-sns-event.json").read_text())
    lambda_handler(event, "")
