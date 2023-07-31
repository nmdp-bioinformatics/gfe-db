import os
import logging
import json
import boto3

logger = logging.getLogger()
logger.setLevel(logging.INFO)

STAGE = os.environ["STAGE"]
APP_NAME = os.environ["APP_NAME"]
AWS_REGION = os.environ["AWS_REGION"]

session = boto3.Session(region_name=AWS_REGION)
ssm = session.client('ssm')

neo4j_backup_document_name = ssm.get_parameter(
    Name=f'/{APP_NAME}/{STAGE}/{AWS_REGION}/Neo4jBackupDocumentName'
)["Parameter"]["Value"]

neo4j_database_instance_id = ssm.get_parameter(
    Name=f'/{APP_NAME}/{STAGE}/{AWS_REGION}/Neo4jDatabaseInstanceId'
)["Parameter"]["Value"]

def lambda_handler(event, context):

    # TODO check if backup is pre or post execution using Parameters field in state machine
    logger.info(json.dumps(event))

    try:
        response = ssm.send_command(
            InstanceIds=[
                neo4j_database_instance_id,
            ],
            DocumentName=neo4j_backup_document_name,
            MaxConcurrency='1',
            CloudWatchOutputConfig={
                'CloudWatchOutputEnabled': True
            })

        if response['ResponseMetadata']['HTTPStatusCode'] != 200:
            logger.error(json.dumps(response, cls=DatetimeEncoder))
            raise Exception("Failed to send command")
        else:
            logger.info(f"Neo4j backup invoked on instance {neo4j_database_instance_id}")

            # TODO poll SSM until command is complete
            # try: poll; except Failed Command: raise Exception
    
    except Exception as err:
        logger.error(err)
        raise err

    return {
        "pre": {
            "DocumentName": response['Command']['DocumentName'],
            "CommandId": response['Command']['CommandId']
        }
    }


# Needed to serialize datetime objects in JSON responses
class DatetimeEncoder(json.JSONEncoder):
    """
    Helps convert datetime objects to pure strings in AWS service API responses. Does not
    convert timezone information.

    Extend `json.JSONEncoder`. 
    """

    def default(self, obj):
        try:
            return super().default(obj)
        except TypeError:
            return str(obj)



if __name__ == "__main__":
    from pathlib import Path

    # event_path = Path(__file__).parent / "pre-execution-event.json"
    event_path = Path(__file__).parent / "post-execution-event.json"

    with open(event_path, "r") as file:
        event = json.load(file)

    lambda_handler(event,"")
