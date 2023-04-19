import logging
from itertools import chain, starmap
from datetime import datetime
import json
import requests
import boto3
from .types import SourceConfig
from .constants import AWS_REGION, GITHUB_PERSONAL_ACCESS_TOKEN

# Logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# boto3 session
session = boto3.Session(region_name=AWS_REGION)
s3 = session.client('s3')

def flatten_json(dictionary, sep='.', skip_fields=[]):
    """Flatten a nested json file. For a list of dictionaries, use this
    inside a for loop before converting to pandas DataFrame."""

    def unpack(parent_key, parent_value):
        """Unpack one level of nesting in json file"""
        # Unpack one level only!!!
        
        if isinstance(parent_value, dict):
            for key, value in parent_value.items():
                temp1 = parent_key + sep + key
                yield temp1, value
        elif isinstance(parent_value, list):
            i = 0 
            for value in parent_value:
                temp2 = parent_key + sep +str(i) 
                i += 1
                yield temp2, value
        else:
            yield parent_key, parent_value    


    # Keep iterating until the termination condition is satisfied
    while True:
        # Keep unpacking the json file until all values are atomic elements (not dictionary or list)
        dictionary = dict(chain.from_iterable(starmap(unpack, dictionary.items())))
        # Terminate condition: not any value in the json file is dictionary or list
        if not any(isinstance(value, dict) for value in dictionary.values()) and \
           not any(isinstance(value, list) for value in dictionary.values()):
            break

    return dictionary


def read_s3_json(bucket, key):
    """Reads config file containing the current state of branches in 
    a GitHub repo"""
    
    try:
        response = s3.get_object(
            Bucket=bucket, 
            Key=key)
        return json.loads(response["Body"].read().decode())
        
    except Exception as err:
        logger.error(f'Failed to read config file to s3://{bucket}/{key}. HTTPStatusCode: {response["ResponseMetadata"]["HTTPStatusCode"]}')
        raise err


def write_s3_json(bucket, key, data):
    """Writes config file containing the current state of branches in 
    a GitHub repo"""
    
    try:
        response = s3.put_object(
            Bucket=bucket, 
            Key=key,
            Body=json.dumps(data).encode())
        
    except Exception as err:
        logger.error(f'Failed to write config file to s3://{bucket}/{key}. HTTPStatusCode: {response["ResponseMetadata"]["HTTPStatusCode"]}')
        raise err


def read_source_config(bucket, key):
    data = read_s3_json(bucket, key)
    return SourceConfig(**data)


def write_source_config(bucket, key, source_config: SourceConfig):
    write_s3_json(bucket, key, source_config.dict())


def get_commits(owner, repo, per_page=100):
    """Return a list of GitHub commits for the specified repository"""

    base_url = 'https://api.github.com'

    # Endpoint
    endpoint = f'/repos/{owner}/{repo}/commits?per_page={per_page}'

    url = base_url + endpoint

    # Headers
    headers = {
        'Authorization': f'token {GITHUB_PERSONAL_ACCESS_TOKEN}',
        'Content-Type': 'application/json',
        'Accept': 'application/vnd.github.v3+json',
        'X-GitHub-Api-Version': '2022-11-28'
    }

    response = requests.get(url, headers=headers)

    return response.json()


def get_commit(owner, repo, commit_sha):
    """Return the commit for the specified repository and commit SHA"""

    base_url = 'https://api.github.com'

    # Endpoint
    endpoint = f'/repos/{owner}/{repo}/commits/{commit_sha}'
    url = base_url + endpoint

    # Headers
    headers = {
        'Authorization': f'token {GITHUB_PERSONAL_ACCESS_TOKEN}',
        'Content-Type': 'application/json',
        'Accept': 'application/vnd.github.v3+json',
        'X-GitHub-Api-Version': '2022-11-28'
    }

    response = requests.get(url, headers=headers)

    return response.json()


def get_file_contents(owner, repo, path):

    base_url = 'https://api.github.com'

    # Endpoint
    endpoint = f'/repos/{owner}/{repo}/contents/{path}'
    url = base_url + endpoint

    # Headers
    headers = {
        'Authorization': f'token {GITHUB_PERSONAL_ACCESS_TOKEN}',
        'Content-Type': 'application/json',
        'Accept': 'application/vnd.github.v3+json',
        'X-GitHub-Api-Version': '2022-11-28'
    }

    response = requests.get(url, headers=headers)

    return response.json()


def get_commits_for_asset(owner, repo, path, since=None):

    base_url = 'https://api.github.com'

    # Endpoint
    endpoint = f'/repos/{owner}/{repo}/commits'
    url = base_url + endpoint

    # Headers
    headers = {
        'Authorization': f'token {GITHUB_PERSONAL_ACCESS_TOKEN}',
        'Content-Type': 'application/json',
        'Accept': 'application/vnd.github.v3+json',
        'X-GitHub-Api-Version': '2022-11-28'
    }

    params = {
        'path': path,
        # validate date is in ISO 8601 format
        'since': since.isoformat() if isinstance(since, datetime) else since
    }

    response = requests.get(url, headers=headers, params=params)

    return response.json()


def get_repo_contents(owner, repo, path, commit_sha=None):

    base_url = 'https://api.github.com'

    # Endpoint
    endpoint = f'/repos/{owner}/{repo}/contents/{path}'
    url = base_url + endpoint

    # Headers
    headers = {
        'Authorization': f'token {GITHUB_PERSONAL_ACCESS_TOKEN}',
        'Content-Type': 'application/json',
        'Accept': 'application/vnd.github.v3+json',
        'X-GitHub-Api-Version': '2022-11-28'
    }

    params = {
        'ref': commit_sha
    }

    response = requests.get(url, headers=headers, params=params)

    return response.json()


def get_repo_asset(owner, repo, path, commit_sha=None):
    """Download a file from a GitHub repository"""
    repo_contents = get_repo_contents(owner, repo, path, commit_sha)

    try:
        res = requests.get(repo_contents['download_url'])
    except KeyError:
        logger.error(json.dumps(repo_contents))
        raise Exception(f'Error downloading {path}')
    
    if res.status_code != 200:
        logger.error(f'Status code {res.status_code} for {path}')
        raise Exception(f'Error downloading {path}')
    
    return res.text


def get_branches(owner, repo):
    """Fetch branches for a GitHub repository"""

    base_url = 'https://api.github.com'

    # Endpoint
    endpoint = f'/repos/{owner}/{repo}/branches'
    url = base_url + endpoint

    # Headers
    headers = {
        'Authorization': f'token {GITHUB_PERSONAL_ACCESS_TOKEN}',
        'Content-Type': 'application/json',
        'Accept': 'application/vnd.github.v3+json',
        'X-GitHub-Api-Version': '2022-11-28'
    }

    response = requests.get(url, headers=headers)
    branches = response.json()

    return branches


def get_branch(owner, repo, branch_name):
    """Fetch branches for a GitHub repository"""

    base_url = 'https://api.github.com'

    # Endpoint
    endpoint = f'/repos/{owner}/{repo}/branches/{branch_name}'
    url = base_url + endpoint

    # Headers
    headers = {
        'Authorization': f'token {GITHUB_PERSONAL_ACCESS_TOKEN}',
        'Content-Type': 'application/json',
        'Accept': 'application/vnd.github.v3+json',
        'X-GitHub-Api-Version': '2022-11-28'
    }

    response = requests.get(url, headers=headers)
    branches = response.json()

    return branches


# Function to fetch pull requests
def get_pull_requests(owner, repo):
    url = f"https://api.github.com/repos/{owner}/{repo}/pulls?state=all"
    
    # Headers
    headers = {
        'Authorization': f'token {GITHUB_PERSONAL_ACCESS_TOKEN}',
        'Content-Type': 'application/json',
        'Accept': 'application/vnd.github.v3+json',
        'X-GitHub-Api-Version': '2022-11-28'
    }
    response = requests.get(url, headers=headers)

    if response.status_code == 200:
        return response.json()
    else:
        print(f"Error: {response.status_code}")
        return []
    

def merge_release_version_with_commit(unique_shas, release_versions):
    # Convert release_versions to a dictionary for easier lookup
    release_versions_dict = {sha: version for sha, version in release_versions}

    # Merge the arrays using a list comprehension
    merged_data = [(release_versions_dict[sha], sha, date) for sha, date in unique_shas]

    # sort by date
    merged_data.sort(key=lambda x: x[2])
    
    return merged_data