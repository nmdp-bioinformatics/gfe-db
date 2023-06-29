import os
import logging
from dataclasses import dataclass
from functools import lru_cache
import re
import json
from typing import Any
import boto3

# Logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


# TODO send to CloudFormation
app_params_mapping = {
    "infra": [
        "VpcID",
        "PublicSubnetID",
        "DataBucketName",
        "DataBucketArn",
        "DataBucketRegionalDomainName",
        "Neo4jDatabaseEndpoint",
        "Neo4jDatabaseEndpointAllocationId",
        "DataPipelineErrorsTopicArn"
    ],
    "pipeline": [
        "GithubSourceRepository",
        "GitHubPersonalAccessToken",
        "ExecutionStateTableName",
        "BuildJobQueueArn",
        "BuildServiceRepositoryName",
        "GfeDbProcessingQueueUrl",
        "GfeDbExecutionResultTopicArn",
        "LoadReleaseActivityArn"
        "UpdatePipelineStateMachineArn",
        "Neo4jLoadQueryDocumentName",
        "DatabaseSyncScriptsDocumentName"
    ],
    "database": [
        "Neo4jCredentialsSecretArn",
        "Neo4jDatabaseSecurityGroupName",
        "Neo4jDatabaseInstanceId",
    ]
}

# TODO send to CloudFormation
app_secrets_mapping = {
    "pipeline": [
        "GitHubPersonalAccessToken"
    ],
    "database": [
        "Neo4jCredentials"
    ]
}

def camel_to_snake(name):
    name = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', name).lower()

# TODO validate that max_size can be total count of all params
@lru_cache
def get_parameter_value(ssm_client, parameter_path: str) -> str:
    logger.info(f"Getting parameter {parameter_path}")
    return ssm_client.get_parameter(
        Name=parameter_path
    )["Parameter"]["Value"]

def get_secret_value(secretsmanager_client, secret_id: str) -> str:
    return secretsmanager_client.get_secret_value(
        SecretId=secret_id
    )["SecretString"]


@dataclass
class JsonModel:
    __slots__ = "__dict__"

    def __init__(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __str__(self) -> str:
        return json.dumps(self.__dict__)

    def __repr__(self) -> str:
        return self.__dict__()


class SessionManager:

    def __init__(self, region_name: str = None) -> None:
        self._region = region_name
        self._session = boto3.Session(region_name=region_name)
        self.clients = JsonModel()  # Cache for clients
        self.resources = JsonModel()  # Cache for resources

    def get_client(self, service_name: str):
        if not hasattr(self.clients, service_name):
            setattr(self.clients, service_name, self._session.client(service_name, region_name=os.environ["AWS_REGION"]))
        return getattr(self.clients, service_name)
    
    def get_resource(self, service_name: str):
        if not hasattr(self.resources, service_name):
            setattr(self.resources, service_name, self._session.resource(service_name, region_name=os.environ["AWS_REGION"]))
        return getattr(self.resources, service_name)
    


class LayerParameterManager(JsonModel):

    def __init__(self, **kwargs) -> None:
        self._params_mapping = kwargs.get("params_mapping")
        self.params_prefix = kwargs.get("params_prefix")
        self._ssm_client = kwargs.get("ssm_client")
    
    def __getattr__(self, name) -> Any:
        if name in self._params_mapping.keys() and name not in self.__dict__.keys():
            value = get_parameter_value(self._ssm_client, f'{self.params_prefix}/{self._params_mapping[name]}')
            setattr(self, name, value)
            return value
        else:
            raise AttributeError(f"Attribute {name} not found")

    def __str__(self) -> str:
        return json.dumps({k: v for k, v in self.__dict__.items() if k != '_ssm_client'})

# TODO
class LayerSecretsManager:
    pass
        

# TODO
class AppLayer:
    def __init__(self, **kwargs) -> None:
        self.layer_name = kwargs.get("layer_name")
        # self.params_prefix = f'/{self.env.APP_NAME}/{self.env.STAGE}/{self.env.AWS_REGION}' # TODO create a named tuple of pydantic class for the path to enforce
        self.params = JsonModel()
        self._init_params()
    
    def _init_params(self):
        setattr(self.params, self.layer_name, LayerParameterManager(
            params_mapping=self.params_mapping[self.layer_name],
            params_prefix=f'{self.params_prefix}', # TODO f'{self.params_prefix}/{layer}',
            ssm_client=self._ssm_client))

class AppConfig:

    def __init__(self, params_mapping: dict, boto3_session = None) -> None:
        self.params_mapping = { k: {camel_to_snake(item): item for item in v} for k, v in app_params_mapping.items() }
        self.env = JsonModel(**{
            "AWS_REGION": os.environ["AWS_REGION"],
            "APP_NAME": os.environ["APP_NAME"],
            "STAGE": os.environ["STAGE"]
        })
        self.session = boto3_session or SessionManager()
        self.params = JsonModel()

        for layer in self.params_mapping:
            setattr(self.params, layer, AppLayer(layer_name=layer, params_mapping=self.params_mapping, boto3_session=self.session))

        

# GITHUB_PERSONAL_ACCESS_TOKEN = secretsmanager.get_secret_value(
#     SecretId=f'/{os.environ["APP_NAME"]}/{os.environ["STAGE"]}/{os.environ["AWS_REGION"]}/GitHubPersonalAccessToken'
# )["SecretString"]


# # TODO put params in Lambda context so that functions that don't need this don't have extra overhead
# # Get SSM Parameters
# github_source_repository = json.loads(ssm.get_parameter(
#     Name=f'/{os.environ["APP_NAME"]}/{os.environ["STAGE"]}/{os.environ["AWS_REGION"]}/GithubSourceRepository'
# )["Parameter"]["Value"])
# GITHUB_REPOSITORY_OWNER = github_source_repository["owner"]
# GITHUB_REPOSITORY_NAME = github_source_repository["name"]

# execution_state_table_name = ssm.get_parameter(
#     Name=f'/{os.environ["APP_NAME"]}/{os.environ["STAGE"]}/{os.environ["AWS_REGION"]}/ExecutionStateTableName'
# )["Parameter"]["Value"]

# # state_machine_arn = ssm.get_parameter(
# #     Name=f'/{os.environ["APP_NAME"]}/{os.environ["STAGE"]}/{os.environ["AWS_REGION"]}/UpdatePipelineArn'
# # )["Parameter"]["Value"]

# data_bucket_name = ssm.get_parameter(
#     Name=f'/{os.environ["APP_NAME"]}/{os.environ["STAGE"]}/{os.environ["AWS_REGION"]}/DataBucketName'
# )["Parameter"]["Value"]

# gfedb_processing_queue_url = ssm.get_parameter(
#     Name=f'/{os.environ["APP_NAME"]}/{os.environ["STAGE"]}/{os.environ["AWS_REGION"]}/GfeDbProcessingQueueUrl'
# )["Parameter"]["Value"]

# TODO move to CloudFormation
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

app = AppConfig(app_params_mapping)
try:
    app.params.pipeline.data_bucket_name
except AttributeError:
    app.params.infra.data_bucket_name
print(app.env)