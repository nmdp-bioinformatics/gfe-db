"""
Lambda with EventBridge event source
* State Machine events trigger the function (State Change, Success, Failure etc.)
* This function may need to have access to the names of states in the state machine so it can update the execution status
  * Map state names to status updates
"""

import os
import logging
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

# Get SSM Parameters
neo4j_database_instance_id = ssm.get_parameter(Name=os.environ["NEO4J_DATABASE_INSTANCE_ID_SSM_PARAM"])["Parameter"]["Value"]
update_pipeline_state_machine_arn = ssm.get_parameter(Name=os.environ["UDPATE_PIPELINE_STATE_MACHINE_ARN_SSM_PARAM"])['Parameter']['Value']

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
    logger.info(json.dumps(event))

    # Invoke state machine

    return

if __name__ == "__main__":
    lambda_handler({}, object)