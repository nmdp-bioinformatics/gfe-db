import os
import logging
import json
import boto3

# set up logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Environment
AWS_REGION = os.environ["AWS_REGION"]

def lambda_handler(event, context):
    logger.info(json.loads(event))

    return