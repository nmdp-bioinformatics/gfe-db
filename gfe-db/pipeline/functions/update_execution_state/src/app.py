import os
import logging
import json
import boto3

# set up logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Environment
AWS_REGION = os.environ["AWS_REGION"]

# Get SSM Parameters
ssm = boto3.client('ssm')
# state_machine_arn '/${AppName}/${Stage}/${AWS::Region}/UpdatePipelineArn'
state_machine_arn = ssm.get_parameter(
    Name=f'/{os.environ["APP_NAME"]}/{os.environ["STAGE"]}/{os.environ["AWS_REGION"]}/UpdatePipelineArn'
)['Parameter']['Value']

# table_name '/${AppName}/${Stage}/${AWS::Region}/ExecutionStateTableName'
table_name = ssm.get_parameter(
    Name=f'/{os.environ["APP_NAME"]}/{os.environ["STAGE"]}/{os.environ["AWS_REGION"]}/ExecutionStateTableName'
)['Parameter']['Value']

# job queue ARN '/${AppName}/${Stage}/${AWS::Region}/BuildJobQueueArn'
job_queue_arn = ssm.get_parameter(
    Name=f'/{os.environ["APP_NAME"]}/{os.environ["STAGE"]}/{os.environ["AWS_REGION"]}/BuildJobQueueArn'
)['Parameter']['Value']

def lambda_handler(event, context):
    logger.info(json.loads(event))

    # Get State Machine ARN

    # Load event

    # Validate Event is from correct State Machine or Batch Job queue

    # Determine event type

    # Update table
    
    return