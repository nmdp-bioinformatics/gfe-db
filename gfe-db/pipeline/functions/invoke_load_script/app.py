import os
import logging
import json
import boto3

logger = logging.getLogger()
logger.setLevel(logging.INFO)

# TODO: Environment variables
data_bucket_name = os.environ["DATA_BUCKET_NAME"]
neo4j_load_query_document_name = os.environ["NEO4J_LOAD_QUERY_DOCUMENT_NAME"]
neo4j_database_instance_id = os.environ["NEO4J_DATABASE_INSTANCE_ID"]

# SSM Document parameters - ssm:GetDocument
ssm = boto3.client('ssm')
response = ssm.get_document(
    Name=neo4j_load_query_document_name)

neo4j_load_query_document_content = json.loads(response["Content"])

# Extract document parameters
neo4j_load_query_document_parameters = neo4j_load_query_document_content["parameters"]
command_line_default = neo4j_load_query_document_parameters["commandLine"]["default"]
source_info_default = neo4j_load_query_document_parameters["sourceInfo"]["default"]

def lambda_handler(event, context):
    """Invoke SSM Run Command for server side loading on Neo4j

    Example using AWS CLI:
    .. code-block:: python
        aws ssm send-command \
        --document-name "dev-gfe-db-database-Neo4jLoadQueryDocument-UgYcOg48yiQB" \
        --document-version "1" \
        --targets '[{"Key":"InstanceIds","Values":["i-0f8ec07e314226283"]}]' \
        --parameters '{"executionTimeout":["3600"],"sourceInfo":["{\"path\":\"https://dev-gfe-db-531868584498-us-east-1.s3.amazonaws.com/config/scripts/load_db.sh\"}"],"sourceType":["S3"],"workingDirectory":["/home/ubuntu"],"commandLine":["bash load_db.sh"]}' \
        --timeout-seconds 600 \
        --max-concurrency "50" \
        --max-errors "0" \
        --cloud-watch-output-config '{"CloudWatchOutputEnabled":true}' \
        --region us-east-1
    
    """

    logger.info(json.dumps(event))
    release = event["RELEASES"]

    # Update params for this execution
    cmd = " ".join([command_line_default, release])

    try:
        response = ssm.send_command(
            InstanceIds=[
                neo4j_database_instance_id,
            ],
            DocumentName=neo4j_load_query_document_name,
            Parameters={
                "commandLine":[cmd],
                "sourceInfo":[json.dumps(source_info_default)]
            },
            MaxConcurrency='1')

        if response['ResponseMetadata']['HTTPStatusCode'] != 200:
            logger.error(json.dumps(response, cls=DatetimeEncoder))
            raise Exception("Failed to send command")
    
    except Exception as err:
        logger.error(err)
        raise err

    return response


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

    path = '/Users/ammon/Documents/00-Projects/nmdp-bioinformatics/02-Repositories/gfe-db/gfe-db/pipeline/functions/invoke_load_script/event.json'

    with open(path, "r") as file:
        event = json.load(file)

    lambda_handler(event,"")
