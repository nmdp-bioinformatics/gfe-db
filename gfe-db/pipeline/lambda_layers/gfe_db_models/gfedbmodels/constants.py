import os
import logging
from typing import Union, Literal
from dataclasses import dataclass
from functools import lru_cache
import re
import json
from typing import Any
import boto3

# Logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


# TODO send to CloudFormation and retrieve all with SSM Parameter "AppConfigLayerMapping"
app_config_layer_mapping = {
    "ssm": {
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
    ]},
    "secretsmanager": {
    "pipeline": [
        "GitHubPersonalAccessToken"
    ],
    "database": [
        "Neo4jCredentials"
    ]
}}

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

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)


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
    
        
class ConfigManager(JsonModel):

    def __init__(self, **kwargs) -> None:
        self._service = kwargs["service"]
        self.map = kwargs["mapping"]
        self._path_prefix = kwargs["path_prefix"]
        self._client = kwargs["client"]
    
    def __getattr__(self, name: str) -> Union[str, dict, list, int, float, bool]:
        """Fetches and caches values from SSM Parameter Store or Secrets Manager. Will not fetch the same value twice.

        Args:
            name (str): Name of the parameter or secret to fetch

        Raises:
            ValueError: Bad value for service type
            AttributeError: Parameter or secret not found

        Returns:
            Union[str, dict, list, int, float, bool]: Value of the parameter or secret
        """
        if name in self.map.keys() and name not in self.__dict__.keys():
            if self._service == "ssm":
                value = get_parameter_value(self._client, f'{self._path_prefix}/{self.map[name]}')
            elif self._service == "secretsmanager":
                value = get_secret_value(self._client, f'{self._path_prefix}/{self.map[name]}')
            else:
                raise ValueError(f"Service {self._service} not supported, must be one of 'ssm' or 'secretsmanager'")
            
            # detect json
            try:
                value = json.loads(value)
            except json.JSONDecodeError:
                pass

            setattr(self, name, value)
            return value
        else:
            raise AttributeError(f"Could not find '{name}' in {self._service} mapping")

    def __str__(self) -> str:
        return json.dumps({k: v for k, v in self.__dict__.items() if k != '_client'})
        

class AppConfig:

    def __init__(self, params_mapping: dict = None, secrets_mapping: dict = None, boto3_session = None, **kwargs) -> None:
        self._build_mappings(params_mapping, secrets_mapping)
        # for name, mapping in {"ssm": params_mapping, "secretsmanager": secrets_mapping}.items():
        #     self.mappings[name] = { k: {camel_to_snake(item): item for item in v} for k, v in mapping.items() }
        self.env = JsonModel(**{
            "AWS_REGION": os.environ["AWS_REGION"],
            "APP_NAME": os.environ["APP_NAME"],
            "STAGE": os.environ["STAGE"]
        })
        self._path_prefix = f'/{self.env.APP_NAME}/{self.env.STAGE}/{self.env.AWS_REGION}' # TODO create a named tuple of pydantic class for the path to enforce
        self.params, self.secrets = JsonModel(), JsonModel() # hard coded model attributes define intended class scope
        self.session = boto3_session or SessionManager()

        for service in ["ssm", "secretsmanager"]:
            setattr(self.session.clients, service, self.session.get_client(service)) 
            self._init_config(service=service)

    def _build_mappings(self, params_mapping: dict, secrets_mapping: dict):
        self.mappings = {}
        for name, mapping in {"ssm": params_mapping, "secretsmanager": secrets_mapping}.items():
            self.mappings[name] = { k: {camel_to_snake(item): item for item in v} for k, v in mapping.items() }

    def _init_config(self, service: str = Union[Literal["ssm"],Literal["secretsmanager"]]):
        for layer in self.mappings[service]:
            setattr(self.params if service == "ssm" else self.secrets, layer, ConfigManager(
                service=service,
                mapping=self.mappings[service][layer],
                path_prefix=f'{self._path_prefix}', # TODO f'{self._path_prefix}/{layer}/Name',
                client=self.session.clients[service]))

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

app = AppConfig(
    params_mapping=app_config_layer_mapping["ssm"],
    secrets_mapping=app_config_layer_mapping["secretsmanager"]
)

# ssm param
try:
    app.params.pipeline.data_bucket_name
except AttributeError:
    app.params.infra.data_bucket_name

# secret
try:
    app.secrets.infra.git_hub_personal_access_token
except AttributeError:
    app.secrets.pipeline.git_hub_personal_access_token
print(app.env)