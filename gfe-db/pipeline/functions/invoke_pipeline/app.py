import os
import logging
import datetime
import copy
import json
import re
import requests
import numpy as np
import boto3

logger = logging.getLogger()
logger.setLevel(logging.INFO)

# TODO: Add database state to JSON that contains all loaded releases
# TODO: make sure that state is being properly updated after each run

GITHUB_PERSONAL_ACCESS_TOKEN = os.environ["GITHUB_PERSONAL_ACCESS_TOKEN"]
GITHUB_REPOSITORY_OWNER = os.environ["GITHUB_REPOSITORY_OWNER"]
GITHUB_REPOSITORY_NAME = os.environ["GITHUB_REPOSITORY_NAME"]
GFE_BUCKET = os.environ["GFE_BUCKET"]
# TODO: add to Makefile
PIPELINE_STATE_PATH = os.environ["PIPELINE_STATE_PATH"]
PIPELINE_PARAMS_PATH = os.environ["PIPELINE_PARAMS_PATH"]
UPDATE_PIPELINE_STATE_MACHINE_ARN = os.environ["UPDATE_PIPELINE_STATE_MACHINE_ARN"]

branches_state_path = f"s3://{GFE_BUCKET}/{PIPELINE_STATE_PATH}"
pipeline_params_path = f"s3://{GFE_BUCKET}/{PIPELINE_PARAMS_PATH}"

s3 = boto3.client('s3')
sfn = boto3.client('stepfunctions')

def lambda_handler(event, context):
    """Checks for new IMGT/HLA releases and triggers the update
    pipeline if any are found"""

    logger.info(json.dumps(event))

    if "releases" in event:
        logging.info(f'Reading parameters from event')
        new_releases, params = parse_event(event)
    else:
        # Load the previous repository state
        logging.info(f'Reading parameters from file')
        new_releases, params = parse_state(branches_state_path, pipeline_params_path)
    
    if new_releases:
        execution_input = []

        for release in new_releases:
            params_input = copy.deepcopy(params)
            params_input["releases"] = release
            params_input = {k.upper():v for k,v in params_input.items()}
            logger.info(f'Running pipeline with these parameters:\n{json.dumps(params_input)}')            
            execution_input.append(params_input)

        payload = {
            "input": execution_input,
            # "error_status": {
            #     "build": False,
            #     "load": False
            # },
        }

        # TODO: include release number in execution identifier
        response = sfn.start_execution(
            stateMachineArn=UPDATE_PIPELINE_STATE_MACHINE_ARN,
            input=json.dumps(payload))

        # Update the config file
        write_config(branches_state_path)

        return {
            # TODO: add timestamp
            "status": response['ResponseMetadata']['HTTPStatusCode'],
            "message": "Pipeline triggered",
            "payload": payload
        }

    else:
        # Update the config file
        write_config(branches_state_path)
        
        return {
            "status": 200,
            "message": "No new releases found"
        }


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


def get_branches(owner, repo):
    """Return a list of GitHub branches for the specified repository"""

    base_url = 'https://api.github.com'

    # Endpoint
    endpoint = f'/repos/{owner}/{repo}/branches?per_page=100'

    url = base_url + endpoint

    # Headers
    headers = {
        'Authorization': f'token {GITHUB_PERSONAL_ACCESS_TOKEN}', 
        'Content-Type': 'application/json',
        'Accept': 'application/vnd.github.v3+json'
    }
    
    response = requests.get(url)
    branches = json.loads(response.content)
    
    return [branch["name"] for branch in branches]


def is_valid_release(branch):
    """Returns True if the branch is a valid release, False if not"""

    # IMGT/HLA release format string
    # Checks for a pattern corresponding to 3 digits followed by one zero, ie., 3460
    release_pattern = r'^\d{3}0$'
    p = re.compile(release_pattern)
    match = p.match(branch)

    if match:
        return True
    else:
        return False


def get_releases(owner, repo):
    """Filters repository branches for those that match the IMGT/HLA release format string"""
    return list(filter(is_valid_release, get_branches(owner, repo)))


def write_config(path):
    """Writes config file containing the current state of branches in 
    a GitHub repo"""
    
    branches_config = {
        "timestamp": str(datetime.datetime.utcnow())[:-7],
        "repository_url": f"https://github.com/{GITHUB_REPOSITORY_OWNER}/{GITHUB_REPOSITORY_NAME}",
        "releases": get_releases(GITHUB_REPOSITORY_OWNER, GITHUB_REPOSITORY_NAME)
    }
    
    try:
        response = s3.put_object(
             Body=json.dumps(branches_config),
             Bucket=GFE_BUCKET,
             Key="/".join(path.split("/")[3:]))
        
        if response['ResponseMetadata']['HTTPStatusCode'] == 200:
            logger.info(f"Wrote config file to {path}")
            return
        else:
            logger.error(f'Failed to write config file to {path}. HTTPStatusCode: {response["ResponseMetadata"]["HTTPStatusCode"]}')
            return
        
    except Exception as err:
        raise err


def read_config(path):
    """Reads config file containing the current state of branches in 
    a GitHub repo"""
    
    try:
        response = s3.get_object(
            Bucket=GFE_BUCKET, 
            Key="/".join(path.split("/")[3:]))
        
        if response['ResponseMetadata']['HTTPStatusCode'] == 200:
            logger.info(f"Read config file from {path}")
            return json.loads(response["Body"].read().decode())
        else:
            logger.error(f'Failed to read config file to {path}. HTTPStatusCode: {response["ResponseMetadata"]["HTTPStatusCode"]}')
            return
        
    except Exception as err:
        raise err


def check_new_releases(previous_state, current_state):
    """Compares the previous repository state as a list of branches, with the current
    repository state. Returns new branches if they are a valid IMGT/HLA release, None 
    if not."""

    # Check if any branches have been added
    new_branches_count = len(current_state) - len(previous_state)
    branches_added = (new_branches_count > 0)

    if branches_added:

        logger.info(f"Found {new_branches_count} new branches: {json.dumps(current_state[-new_branches_count:])}")

        # Get the new branches
        new_releases = sorted([int(release) for release in list(set(current_state).difference(previous_state))])
        previous_state_last_release = [int(previous_state[-1])]

        # Check that the last release of the previous state and the newly added releases all differ by 10
        elementwise_difference = list(set(np.diff([release for release in previous_state_last_release + new_releases])))
        new_branches_are_valid_releases = (len(elementwise_difference) == 1 and elementwise_difference[0] == 10)

        if new_branches_are_valid_releases:
            
            return [str(release) for release in new_releases]
    else:
        logger.info("No new branches detected")
        
        return

def check_current_executions(state_machine_arn):
    
    response = sfn.list_executions(
        stateMachineArn=state_machine_arn,
        statusFilter='RUNNING')

    # Extract executions
    executions_arns = [execution['executionArn'] for execution in response['executions']]
    
    releases_processing = []

    for executions_arn in executions_arns:

        response = sfn.describe_execution(
            executionArn=executions_arn)

        releases_processing = releases_processing + [params["releases"] for params in json.loads(response['input'])]

    return releases_processing
    

def parse_event(event):
    """Restructures the event and returns pipeline execution parameters"""

    new_releases = str(event["releases"]).replace(" ", "").split(",")
    params = {k:v for k,v in event.items() if k != "releases"}

    return new_releases, params


def parse_state(state_path, params_path):
    """"""
    previous_state = read_config(state_path)['releases']

    # Get the current repository state
    current_state = get_releases(GITHUB_REPOSITORY_OWNER, GITHUB_REPOSITORY_NAME)

    # Log repository states
    logger.info(f'Previous repository state:\n{json.dumps(previous_state)}')
    logger.info(f'Current repository state:\n{json.dumps(current_state)}')

    # Compare the previous state with the current state
    new_releases = check_new_releases(previous_state, current_state)

    # TODO: Describe current executions and make sure the release is not already being built
    # releases_processing = set(check_current_executions(UPDATE_PIPELINE_STATE_MACHINE_ARN))

    # Load the current pipeline params
    params = read_config(params_path)

    return new_releases, params


if __name__ == "__main__":
    import os
    from pathlib import Path

    def read_json(path):
        with open(path, 'r') as file:
            event = json.load(file)
        return event

    path = Path(__file__)
    path = str(path.parent.parent / "events" / "event.json")
    
    event = read_json(path)

    lambda_handler("","")
