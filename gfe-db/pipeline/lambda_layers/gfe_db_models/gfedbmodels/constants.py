import os
import logging
import json
import boto3

# Logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Environment variables
AWS_REGION = os.environ["AWS_REGION"] 

# AWS clients
session = boto3.Session(region_name=AWS_REGION)
ssm = session.client("ssm")
dynamodb = session.resource("dynamodb", region_name=AWS_REGION)
sqs = session.client("sqs", region_name=AWS_REGION)
secretsmanager = session.client("secretsmanager", region_name=AWS_REGION)

GITHUB_PERSONAL_ACCESS_TOKEN = secretsmanager.get_secret_value(
    SecretId=f'/{os.environ["APP_NAME"]}/{os.environ["STAGE"]}/{os.environ["AWS_REGION"]}/GitHubPersonalAccessToken'
)["SecretString"]

# Get SSM Parameters
github_source_repository = json.loads(ssm.get_parameter(
    Name=f'/{os.environ["APP_NAME"]}/{os.environ["STAGE"]}/{os.environ["AWS_REGION"]}/GithubSourceRepository'
)["Parameter"]["Value"])
GITHUB_REPOSITORY_OWNER = github_source_repository["owner"]
GITHUB_REPOSITORY_NAME = github_source_repository["name"]

execution_state_table_name = ssm.get_parameter(
    Name=f'/{os.environ["APP_NAME"]}/{os.environ["STAGE"]}/{os.environ["AWS_REGION"]}/ExecutionStateTableName'
)["Parameter"]["Value"]

# state_machine_arn = ssm.get_parameter(
#     Name=f'/{os.environ["APP_NAME"]}/{os.environ["STAGE"]}/{os.environ["AWS_REGION"]}/UpdatePipelineArn'
# )["Parameter"]["Value"]

data_bucket_name = ssm.get_parameter(
    Name=f'/{os.environ["APP_NAME"]}/{os.environ["STAGE"]}/{os.environ["AWS_REGION"]}/DataBucketName'
)["Parameter"]["Value"]

gfedb_processing_queue_url = ssm.get_parameter(
    Name=f'/{os.environ["APP_NAME"]}/{os.environ["STAGE"]}/{os.environ["AWS_REGION"]}/GfeDbProcessingQueueUrl'
)["Parameter"]["Value"]

# list of fields to be used in the execution state table
execution_state_table_fields = [
    "commit.sha",
    "execution.version",
    "commit.date_utc",
    "commit.html_url",
    "commit.message",
    "execution.date_utc",
    "execution.status",
    "execution.input_parameters.align",
    "execution.input_parameters.kir",
    "execution.input_parameters.limit",
    "execution.input_parameters.mem_profile",
    "repository.name",
    "repository.owner",
    "repository.url",
    "created_utc",
    "updated_utc"
]