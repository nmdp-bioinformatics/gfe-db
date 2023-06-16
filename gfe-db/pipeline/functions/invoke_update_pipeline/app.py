import os
import logging
from datetime import datetime
import json
import boto3

# set up logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Environment
AWS_REGION = os.environ["AWS_REGION"]

# Boto3 Clients
session = boto3.Session(region_name=AWS_REGION)
ssm = session.client('ssm', region_name=AWS_REGION)
ec2 = session.client('ec2', region_name=AWS_REGION)
states = session.client('stepfunctions', region_name=AWS_REGION)
sqs = session.client('sqs', region_name=AWS_REGION)

# Get SSM Parameters
neo4j_database_instance_id = ssm.get_parameter(Name=os.environ["NEO4J_DATABASE_INSTANCE_ID_SSM_PARAM"])["Parameter"]["Value"]
update_pipeline_state_machine_arn = ssm.get_parameter(Name=os.environ["UDPATE_PIPELINE_STATE_MACHINE_ARN_SSM_PARAM"])['Parameter']['Value']
gfe_db_processing_queue_url = ssm.get_parameter(Name=os.environ["GFE_DB_PROCESSING_QUEUE_URL_SSM_PARAM"])['Parameter']['Value']

# Check that database is running, abort if not
try:
    response = ec2.describe_instance_status(InstanceIds=[neo4j_database_instance_id])
    if response['InstanceStatuses'][0]['InstanceState']['Name'] != 'running':
        raise Exception(f'Instance {neo4j_database_instance_id} is not running, aborting...')
    else:
        logger.info(f'Instance {neo4j_database_instance_id} is running')
except Exception as e:
    raise e


def lambda_handler(event, context):
    
    errors = 0
    execution_arns = []
    for record in event['Records']:

        try:
            message = json.loads(record['body'])
            logger.info(f"Received message for version {message['version']} and commit {message['commit_sha']}")
            response = states.start_execution(
                stateMachineArn=update_pipeline_state_machine_arn,
                name=generate_execution_id(message),
                input=json.dumps(message)
            )

            execution_arns.append(response['executionArn'])
            
            try:
                response = sqs.delete_message(
                    QueueUrl=gfe_db_processing_queue_url,
                    ReceiptHandle=record['receiptHandle']
                )
                logger.info(f"Message deleted from queue")
            except Exception as e:
                logger.error(f"Error deleting message from queue: {e}")

        except Exception as e:
            import traceback
            msg = f'Error processing commit {message["commit_sha"]}: {e}\n{traceback.format_exc()}'
            logger.error(msg)
            errors += 1
            continue

    return_msg = f'{len(event["Records"])-errors} of {len(event["Records"])} messages processed successfully, {errors} error(s)'
    if errors > 0:
        logger.error(json.dumps({
            'message': return_msg,
            'execution_arns': execution_arns
        }))
        logger.error(json.dumps(event))
        raise Exception(return_msg)

    return {
        'statusCode': 200,
        'body': json.dumps({
            'message': return_msg,
            'execution_arns': execution_arns
        })
    }


def generate_execution_id(message: dict) -> str:
    """Generate an execution id for the state machine execution with format:
    {version}_{commit_sha}_{YYYYMMDD_HHMMSS}

    Args:
        message (dict): Message from SQS queue

    Returns:
        str: Execution id
    """
    return "_".join([str(message['version']), message['commit_sha'], datetime.utcnow().strftime("%y%m%d_%H%M%S")])


if __name__ == "__main__":
    from pathlib import Path

    with open(Path(__file__).parent / "sqs-event.json", "r") as f:
        event = json.load(f)
    
    lambda_handler(event, "")