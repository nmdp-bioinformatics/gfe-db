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
    pipeline,
    database
)

logger = logging.getLogger()
logger.setLevel(logging.INFO)

neo4j_load_query_document_name = pipeline.params.Neo4jLoadQueryDocumentName
neo4j_database_instance_id = database.params.Neo4jDatabaseInstanceId

# Get SSM Document Neo4jLoadQuery
ssm = session.clients["ssm"]
response = ssm.get_document(Name=neo4j_load_query_document_name)
neo4j_load_query_document_content = json.loads(response["Content"])

# Extract document parameters
neo4j_load_query_document_parameters = neo4j_load_query_document_content["parameters"]
source_info_default = neo4j_load_query_document_parameters["sourceInfo"]["default"]


def lambda_handler(event, context):
    """Invoke SSM Run Command for server side loading on Neo4j

    Example using AWS CLI:
    .. code-block:: python
        aws ssm send-command \
        --document-name "dev-gfe-db-database-Neo4jLoadQueryDocument-UgYcOg48yiQB" \
        --document-version "1" \
        --targets '[{"Key":"InstanceIds","Values":["i-0f8ec07e314226283"]}]' \
        --parameters '{"executionTimeout":["3600"],"sourceInfo":["{\"path\":\"https://<data bucket name>.s3.amazonaws.com/config/database/scripts/load_db.sh\"}"],"sourceType":["S3"],"workingDirectory":["/home/ec2-user"],"LoadEvent":["{\"key\":\"value\"}"]}' \
        --timeout-seconds 600 \
        --max-concurrency "50" \
        --max-errors "0" \
        --cloud-watch-output-config '{"CloudWatchOutputEnabled":true}' \
        --region us-east-1
    
    """

    logger.info(json.dumps(event))

    # TODO BOOKMARK 5/31/23: Check if Neo4jLoadQueryDocument is already running, if it is exit 0 (makes service idempotent)
    # Note: Neo4jLoadQueryDocument only needs to be triggered once and it will fetch the next release until there are no more left

    try:
        response = ssm.send_command(
            InstanceIds=[
                neo4j_database_instance_id,
            ],
            DocumentName=neo4j_load_query_document_name,
            Parameters={
                "sourceType": ["S3"],
                "sourceInfo": [json.dumps(source_info_default)],
                "workingDirectory": ["/home/ec2-user"],
                "executionTimeout": ["28800"],
                "LoadEvent": [json.dumps(event)],
            },
            MaxConcurrency="1",
            CloudWatchOutputConfig={"CloudWatchOutputEnabled": True},
        )

        if response["ResponseMetadata"]["HTTPStatusCode"] != 200:
            logger.error(json.dumps(response, cls=DatetimeEncoder))
            message = f"Failed to send command to instance {neo4j_database_instance_id}"
            raise Exception("Failed to send command")
        else:
            message = f"Command invoked on instance {neo4j_database_instance_id}"
            logger.info(message)
    
            return {
                "message": message,
                "sqs": event["sqs"],
                "ssm": {
                    "CommandId": response["Command"]["CommandId"],
                    "InstanceId": neo4j_database_instance_id,
                }
            }
        
    except Exception as err:
        logger.error(err)
        raise err



# Serializes datetime objects in JSON responses
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

    event_path = Path(__file__).parent / "event.json"

    with open(event_path, "r") as file:
        event = json.load(file)

    lambda_handler(event, "")
