import os
import logging
from dotenv import load_dotenv, find_dotenv
load_dotenv(find_dotenv());

# Logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Environment variables
AWS_REGION = os.environ["AWS_REGION"] 
GITHUB_PERSONAL_ACCESS_TOKEN = os.environ["GITHUB_PERSONAL_ACCESS_TOKEN"]