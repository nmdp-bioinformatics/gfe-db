import os
import logging
import datetime
import json
import requests
from dotenv import load_dotenv, find_dotenv

load_dotenv(find_dotenv());

logger = logging.getLogger()
logger.setLevel(logging.INFO)

GITHUB_PERSONAL_ACCESS_TOKEN = os.environ["GITHUB_PERSONAL_ACCESS_TOKEN"]
GITHUB_REPOSITORY_OWNER = os.environ["GITHUB_REPOSITORY_OWNER"]
GITHUB_REPOSITORY_NAME = os.environ["GITHUB_REPOSITORY_NAME"]

# Time in seconds between scheduled invocations
LAMBDA_INVOCATION_INTERVAL_SECONDS = int(os.environ["LAMBDA_INVOCATION_INTERVAL_SECONDS"])

base_url = 'https://api.github.com'

headers = {
    'Authorization': f'token {GITHUB_PERSONAL_ACCESS_TOKEN}',
    'Content-Type': 'application/json',
    'Accept': 'application/vnd.github.v3+json'
}

def lambda_handler(event, context):

    # Get latest release version
    endpoint = f'/repos/{GITHUB_REPOSITORY_OWNER}/{GITHUB_REPOSITORY_NAME}/branches?per_page=100'
    url = "".join([base_url, endpoint])
    response = requests.get(url)
    latest_release_version = json.loads(response.content)[-2]['name']

    logger.info(f"latest_release_version: {latest_release_version}")

    # Get most recent commit sha and date
    endpoint = f'/repos/{GITHUB_REPOSITORY_OWNER}/{GITHUB_REPOSITORY_NAME}/commits?per_page=100'
    url = "".join([base_url, endpoint])
    response = requests.get(url)
    commits = json.loads(response.content)

    logger.info(f"commits: {commits}")

    most_recent_commit_sha = commits[0]['sha']
    most_recent_commit_date = commits[0]['commit']['author']['date']
    most_recent_commit_date = datetime.datetime.strptime(most_recent_commit_date, "%Y-%m-%dT%H:%M:%SZ")
    utc_now = datetime.datetime.utcnow()
    time_delta = utc_now - most_recent_commit_date

    logger.info(f"commits: {commits}")

    # If the time passed since the last commit is less than 3600 seconds, trigger the pipeline
    run_pipeline = time_delta.total_seconds() < LAMBDA_INVOCATION_INTERVAL_SECONDS

    # List files from most recent commit
    endpoint = f'/repos/{GITHUB_REPOSITORY_OWNER}/{GITHUB_REPOSITORY_NAME}/contents?ref={most_recent_commit_sha}'
    url = "".join([base_url, endpoint])
    response = requests.get(url)
    contents = json.loads(response.content)
    hla_dat = list(filter(lambda item: item['name'] == 'hla.dat', contents))

    # if hla.dat is in the most recent commit, check if the commit occurred after the last invocation
    if run_pipeline:
        print("Running pipeline...")

    return


# if __name__ == "__main__":
#     lambda_handler("","")
