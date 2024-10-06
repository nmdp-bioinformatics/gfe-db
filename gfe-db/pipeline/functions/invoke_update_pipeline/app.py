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
    pipeline,
    database
)

# set up logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Environment
AWS_REGION = os.environ["AWS_REGION"]

# Boto3 Clients
ec2 = session.client("ec2", region_name=AWS_REGION)
states = session.client("stepfunctions", region_name=AWS_REGION)
sqs = session.client("sqs", region_name=AWS_REGION)

# Get SSM Parameters
neo4j_database_instance_id = database.params.Neo4jDatabaseInstanceId
update_pipeline_state_machine_arn = pipeline.params.UpdatePipelineStateMachineArn
gfe_db_processing_queue_url = pipeline.params.GfeDbProcessingQueueUrl

# Check that database is running, abort if not
# TODO State table should be updated to FAILED if database is not running, however this block needs to run inside the handler
response = ec2.describe_instance_status(InstanceIds=[neo4j_database_instance_id])
if len(response["InstanceStatuses"]) > 0:
    if response["InstanceStatuses"][0]["InstanceState"]["Name"] != "running":
        raise Exception(
            f"Instance {neo4j_database_instance_id} is not running, aborting..."
        )
    else:
        logger.info(f"Instance {neo4j_database_instance_id} is running")
else:
    raise Exception(f"Instance {neo4j_database_instance_id} not found, aborting...")


def lambda_handler(event, context):
    errors = 0
    execution_arns = []
    for record in event["Records"]:
        try:
            message = json.loads(record["body"])

            # Include receipt handle in message to return to queue if step function fails
            message["receipt_handle"] = record["receiptHandle"]
            
            logger.info(
                f"Received message for version {message['version']} and commit {message['commit_sha']}"
            )

            payload = {
                "input": message
            }
            response = states.start_execution(
                stateMachineArn=update_pipeline_state_machine_arn,
                name=payload["input"]["id"], # {version}_{commit_sha}_{YYYYMMDD_HHMMSS}
                input=json.dumps(payload),
            )

            execution_arns.append(response["executionArn"])

            # try:
            #     response = sqs.delete_message(
            #         QueueUrl=gfe_db_processing_queue_url,
            #         ReceiptHandle=record["receiptHandle"],
            #     )
            #     logger.info(f"Message deleted from queue")
            # except Exception as e:
            #     logger.error(f"Error deleting message from queue: {e}")

        except Exception as e:
            import traceback
            message = f'Error processing commit {message["commit_sha"]}: {e}\n{traceback.format_exc()}'
            logger.error(message)
            errors += 1
            continue

    return_msg = f'{len(event["Records"])-errors} of {len(event["Records"])} messages processed successfully, {errors} error(s)'
    if errors > 0:
        message = json.dumps({"message": return_msg, "execution_arns": execution_arns})
        logger.error(message)
        logger.error(json.dumps(event))
        raise Exception(return_msg)
    else:
        message = json.dumps({"message": return_msg, "execution_arns": execution_arns})
        return {
            "statusCode": 200,
            "body": json.dumps({"message": return_msg, "execution_arns": execution_arns}),
        }

# # TODO move to CheckSourceUpdate
# def generate_execution_id(message: dict) -> str:
#     """Generate an execution id for the state machine execution with format:
#     {version}_{commit_sha}_{YYYYMMDD_HHMMSS}

#     Args:
#         message (dict): Message from SQS queue

#     Returns:
#         str: Execution id
#     """
#     return "_".join(
#         [
#             str(message["version"]),
#             message["commit_sha"],
#             datetime.utcnow().strftime("%y%m%d_%H%M%S"),
#         ]
#     )


if __name__ == "__main__":
    from pathlib import Path

    with open(Path(__file__).parent / "sqs-event.json", "r") as f:
        event = json.load(f)

    lambda_handler(event, "")
