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

AWS_REGION = os.environ["AWS_REGION"]
APP_NAME = os.environ["APP_NAME"]
STAGE = os.environ["STAGE"]

# TODO send to CloudFormation and retrieve all with SSM Parameter "AppConfigLayerMapping"
infra_mapping = {
    "ssm": [
            f"/{APP_NAME}/{STAGE}/{AWS_REGION}/VpcID",
            f"/{APP_NAME}/{STAGE}/{AWS_REGION}/PublicSubnetID",
            f"/{APP_NAME}/{STAGE}/{AWS_REGION}/DataBucketName",
            f"/{APP_NAME}/{STAGE}/{AWS_REGION}/DataBucketArn",
            f"/{APP_NAME}/{STAGE}/{AWS_REGION}/DataBucketRegionalDomainName",
            f"/{APP_NAME}/{STAGE}/{AWS_REGION}/Neo4jDatabaseEndpoint",
            f"/{APP_NAME}/{STAGE}/{AWS_REGION}/Neo4jDatabaseEndpointAllocationId",
            f"/{APP_NAME}/{STAGE}/{AWS_REGION}/DataPipelineErrorsTopicArn",
        ]
}

pipeline_mapping = {
    "ssm": [
        f"/{APP_NAME}/{STAGE}/{AWS_REGION}/GithubSourceRepository",
        f"/{APP_NAME}/{STAGE}/{AWS_REGION}/GitHubPersonalAccessToken",
        f"/{APP_NAME}/{STAGE}/{AWS_REGION}/ExecutionStateTableName",
        f"/{APP_NAME}/{STAGE}/{AWS_REGION}/BuildJobQueueArn",
        f"/{APP_NAME}/{STAGE}/{AWS_REGION}/BuildServiceRepositoryName",
        f"/{APP_NAME}/{STAGE}/{AWS_REGION}/GfeDbProcessingQueueUrl",
        f"/{APP_NAME}/{STAGE}/{AWS_REGION}/GfeDbExecutionResultTopicArn",
        f"/{APP_NAME}/{STAGE}/{AWS_REGION}/LoadReleaseActivityArn",
        f"/{APP_NAME}/{STAGE}/{AWS_REGION}/UpdatePipelineStateMachineArn",
        f"/{APP_NAME}/{STAGE}/{AWS_REGION}/Neo4jLoadQueryDocumentName",
        f"/{APP_NAME}/{STAGE}/{AWS_REGION}/DatabaseSyncScriptsDocumentName",
    ],
    "secretsmanager": [
        f"/{APP_NAME}/{STAGE}/{AWS_REGION}/GitHubPersonalAccessToken"
    ],
}

database_mapping = {
    "ssm": [
        f"/{APP_NAME}/{STAGE}/{AWS_REGION}/Neo4jCredentialsSecretArn",
        f"/{APP_NAME}/{STAGE}/{AWS_REGION}/Neo4jDatabaseSecurityGroupName",
        f"/{APP_NAME}/{STAGE}/{AWS_REGION}/Neo4jDatabaseInstanceId",
    ],
    "secretsmanager":  [
        f"/{APP_NAME}/{STAGE}/{AWS_REGION}/Neo4jCredentials"
    ]
}


def camel_to_snake(name):
    name = re.sub("(.)([A-Z][a-z]+)", r"\1_\2", name)
    return re.sub("([a-z0-9])([A-Z])", r"\1_\2", name).lower()


@lru_cache
def get_parameter_value(ssm_client, parameter_path: str) -> str:
    logger.info(f"Getting parameter {parameter_path}")
    return ssm_client.get_parameter(Name=parameter_path)["Parameter"]["Value"]


@lru_cache
def get_secret_value(secretsmanager_client, secret_id: str) -> str:
    return secretsmanager_client.get_secret_value(SecretId=secret_id)["SecretString"]


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
    def __init__(self, boto3_session=None, region_name: str = None) -> None:
        self.session = boto3_session or boto3.Session(region_name=region_name)
        self.clients = JsonModel()  # Cache for clients
        self.resources = JsonModel()  # Cache for resources

    def get_client(self, service_name: str):
        if not hasattr(self.clients, service_name):
            setattr(self.clients, service_name, self.session.client(service_name))
        return getattr(self.clients, service_name)

    def get_resource(self, service_name: str):
        if not hasattr(self.resources, service_name):
            setattr(self.resources, service_name, self.session.resource(service_name))
        return getattr(self.resources, service_name)

    def __getattr__(self, name):
        return getattr(self.session, name)


class ConfigManager(JsonModel):
    def __init__(self, **kwargs) -> None:
        self._service = kwargs["service"]
        self.map = kwargs["mapping"]
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
                value = get_parameter_value(
                    self._client, self.map[name]
                )
            elif self._service == "secretsmanager":
                value = get_secret_value(
                    self._client, self.map[name]
                )
            else:
                raise ValueError(
                    f"Service {self._service} not supported, must be one of 'ssm' or 'secretsmanager'"
                )

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
        return json.dumps({k: v for k, v in self.__dict__.items() if k != "_client"})


class AppConfig:
    """Class to manage the configuration of the application using SSM Parameter Store and Secrets Manager.

    Usage:
        >>> config = AppConfig(ssm_mapping=ssm_mapping, secrets_mapping=secrets_mapping)
        >>> # Show application environment variables
        >>> config.env
        >>> # Fetch a parameter
        >>> config.pipeline.github_source_repository
        >>> # Fetch a secret
        >>> config.database.neo4j_credentials

    Args:
        ssm_mapping (dict, optional): Mapping of parameter names to their path in SSM Parameter Store. Defaults to None.
        secrets_mapping (dict, optional): Mapping of secret names to their path in Secrets Manager. Defaults to None.
        boto3_session ([type], optional): Boto3 session to use. Defaults to None.
    """

    def __init__(
        self,
        ssm_paths: dict = None,
        secrets_paths: dict = None,
        path_separator: str = "/",
        boto3_session=None,
        region_name: str = None,
    ) -> None:
        self.services, self.mappings = self._build_mappings(
            ssm_paths=ssm_paths, 
            secrets_paths=secrets_paths,
            path_separator=path_separator,
        )
        self.session = SessionManager(boto3_session, region_name)

        for service in self.services:
            setattr(self.session.clients, service, self.session.get_client(service))
            setattr(
                self,
                "params" if service == "ssm" else "secrets",
                ConfigManager(
                    service=service,
                    mapping=self.mappings[service],
                    client=self.session.clients[service],
                )
            )

    @staticmethod
    def _build_mappings(**kwargs):
        collections = {
            k: v
            for k, v in {
                "ssm": kwargs.get("ssm_paths"),
                "secretsmanager": kwargs.get("secrets_paths"),
            }.items()
            if v is not None
        }
        services = list(collections.keys())
        mappings = {}
        for name, mapping in collections.items():
            mappings[name] = \
            {
                camel_to_snake(item.split(kwargs.get("path_separator"))[-1]): item for item in mapping
            }
        return services, mappings


# TODO move to CloudFormation
# list of fields to be used in the execution state table
# execution_state_table_fields = [
#     "commit.sha",
#     "execution.version",
#     "commit.date_utc",
#     "commit.html_url",
#     "commit.message",
#     "execution.date_utc",
#     "execution.status",
#     "execution.input_parameters.align",
#     "execution.input_parameters.kir",
#     "execution.input_parameters.limit",
#     "execution.input_parameters.mem_profile",
#     "repository.name",
#     "repository.owner",
#     "repository.url",
#     "created_utc",
#     "updated_utc",
# ]

# session = boto3.Session(region_name="us-east-1")

pipeline = AppConfig(
    ssm_paths=pipeline_mapping["ssm"],
    secrets_paths=pipeline_mapping["secretsmanager"],
    # boto3_session=session
)

print(pipeline.params.map)
# print(pipeline.params.github_source_repository)

# TODO implement AppConfigLayerMapping and test
# Notes
# Maintains least-privileged permissions to parameters and secrets
# minimize configuration
# immediate access to available paramters and secrets