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