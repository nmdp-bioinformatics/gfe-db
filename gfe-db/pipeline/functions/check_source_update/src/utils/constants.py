import os
import logging

# Logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Environment variables
AWS_REGION = os.environ["AWS_REGION"] 
GITHUB_PERSONAL_ACCESS_TOKEN = os.environ["GITHUB_PERSONAL_ACCESS_TOKEN"]
GITHUB_REPOSITORY_OWNER = os.environ["GITHUB_REPOSITORY_OWNER"]
GITHUB_REPOSITORY_NAME = os.environ["GITHUB_REPOSITORY_NAME"]

# list of fields to be used in the execution state table
execution_state_table_fields = [
    "commit.sha",
    "execution.version",
    "commit.date_utc",
    "commit.html_url",
    "commit.message",
    "execution.date_utc",
    "execution.input_parameters",
    "execution.status",
    "repository.default_input_parameters.align",
    "repository.default_input_parameters.kir",
    "repository.default_input_parameters.limit",
    "repository.default_input_parameters.mem_profile",
    "repository.name",
    "repository.owner",
    "repository.url"
]