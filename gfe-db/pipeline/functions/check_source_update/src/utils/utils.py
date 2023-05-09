import os
import logging
from typing import List, Dict
from pathlib import Path
from itertools import chain, starmap
from datetime import datetime
import json
import pickle
import re
import requests
import multiprocessing
from concurrent.futures import ThreadPoolExecutor, as_completed
import boto3
from botocore.exceptions import ClientError
from .types import (
    RepositoryConfig,
    TargetMetadataConfig,
    Commit,
    ExecutionStateItem,
    ExecutionDetailsConfig
)
from .constants import (
    AWS_REGION,
    GITHUB_PERSONAL_ACCESS_TOKEN,
    GITHUB_REPOSITORY_OWNER,
    GITHUB_REPOSITORY_NAME,
)

# Logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# boto3 session
session = boto3.Session(region_name=AWS_REGION)
s3 = session.client("s3")

cache_dir = Path(__file__).parent / "_cache"


def save_json_to_cache(data, var_name):
    """Saves data to cache directory"""
    if not cache_dir.exists():
        cache_dir.mkdir()
    with open(cache_dir / var_name, "w") as f:
        json.dump(data, f, indent=4)


def save_pickle_to_cache(data, var_name):
    """Saves data to cache directory"""
    if not cache_dir.exists():
        cache_dir.mkdir()
    with open(cache_dir / var_name, "wb") as f:
        pickle.dump(data, f)


def load_json_from_cache(var_name):
    """Loads data from cache directory"""
    with open(cache_dir / var_name, "r") as f:
        data = json.load(f)
    return data


def load_pickle_from_cache(var_name):
    """Loads data from cache directory"""
    with open(cache_dir / var_name, "rb") as f:
        data = pickle.load(f)
    return data


# implement a @cache_json decorator to cache the results of the function in a file or load from cache if it exists
def cache_json(func):
    """Decorator to cache function results"""

    def wrapper(*args, **kwargs):
        var_name = func.__name__
        if (cache_dir / var_name).exists():
            logger.info(f"Loading {var_name} from cache")
            return load_json_from_cache(var_name)
        else:
            logger.info(f"Saving {var_name} to cache")
            data = func(*args, **kwargs)
            save_json_to_cache(data, var_name)
            return data

    return wrapper


# rewrite the cache_json decorator to work for pickle files
def cache_pickle(func):
    """Decorator to cache function results"""

    def wrapper(*args, **kwargs):
        var_name = func.__name__
        if (cache_dir / var_name).exists():
            logger.info(f"Loading {var_name} from cache")
            return load_pickle_from_cache(var_name)
        else:
            logger.info(f"Saving {var_name} to cache")
            data = func(*args, **kwargs)
            save_pickle_to_cache(data, var_name)
            return data

    return wrapper


def flatten_json(dictionary, sep=".", skip_fields=[]):
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
                temp2 = parent_key + sep + str(i)
                i += 1
                yield temp2, value
        else:
            yield parent_key, parent_value

    # Keep iterating until the termination condition is satisfied
    while True:
        # Keep unpacking the json file until all values are atomic elements (not dictionary or list)
        dictionary = dict(chain.from_iterable(starmap(unpack, dictionary.items())))
        # Terminate condition: not any value in the json file is dictionary or list
        if not any(
            isinstance(value, dict) for value in dictionary.values()
        ) and not any(isinstance(value, list) for value in dictionary.values()):
            break

    if len(skip_fields) > 0:
        return {k: v for k, v in dictionary.items() if k not in skip_fields}
    else:
        return dictionary


def read_s3_json(bucket, key):
    """Reads config file containing the current state of branches in
    a GitHub repo"""

    try:
        response = s3.get_object(Bucket=bucket, Key=key)
        return json.loads(response["Body"].read().decode())

    except ClientError as err:
        logger.error(f"Failed to read config file to s3://{bucket}/{key}")
        raise err


def write_s3_json(bucket, key, data):
    """Writes config file containing the current state of branches in
    a GitHub repo"""

    try:
        response = s3.put_object(Bucket=bucket, Key=key, Body=json.dumps(data).encode())

    except Exception as err:
        logger.error(
            f'Failed to write config file to s3://{bucket}/{key}. HTTPStatusCode: {response["ResponseMetadata"]["HTTPStatusCode"]}'
        )
        raise err


# def read_source_config(bucket, key):
#     data = read_s3_json(bucket, key)
#     return SourceConfig(**data)


# def write_source_config(bucket, key, source_config: SourceConfig):
#     write_s3_json(bucket, key, source_config.dict())


def list_commits(owner, repo, **kwargs):
    """Return a list of GitHub commits for the specified repository"""

    base_url = "https://api.github.com"

    # Endpoint
    endpoint = f"/repos/{owner}/{repo}/commits"

    url = base_url + endpoint

    params = {
        "per_page": kwargs.get("per_page"),
        "page": kwargs.get("page"),
    }

    # Headers
    headers = {
        "Authorization": f"token {GITHUB_PERSONAL_ACCESS_TOKEN}",
        "Content-Type": "application/json",
        "Accept": "application/vnd.github.v3+json",
        "X-GitHub-Api-Version": "2022-11-28",
    }

    response = requests.get(url, headers=headers, params=params)

    return response.json()


# @cache_json
def paginate_commits(owner, repo, start_page=1, per_page=100, **kwargs):
    page = start_page
    commits = []
    while True:
        response = list_commits(owner, repo, page=page, per_page=per_page, **kwargs)
        if len(response) == 0:
            break
        logger.debug(f"Page {page}: {len(response)} commits")
        commits.extend(response)
        page += 1

    if len(commits) == 0:
        raise ValueError("No commits found")

    return commits


def get_commit(owner, repo, commit_sha):
    """Return the commit for the specified repository and commit SHA"""

    base_url = "https://api.github.com"

    # Endpoint
    endpoint = f"/repos/{owner}/{repo}/commits/{commit_sha}"
    url = base_url + endpoint

    # Headers
    headers = {
        "Authorization": f"token {GITHUB_PERSONAL_ACCESS_TOKEN}",
        "Content-Type": "application/json",
        "Accept": "application/vnd.github.v3+json",
        "X-GitHub-Api-Version": "2022-11-28",
    }

    response = requests.get(url, headers=headers)

    return response.json()


def get_file_contents(owner, repo, path):
    base_url = "https://api.github.com"

    # Endpoint
    endpoint = f"/repos/{owner}/{repo}/contents/{path}"
    url = base_url + endpoint

    # Headers
    headers = {
        "Authorization": f"token {GITHUB_PERSONAL_ACCESS_TOKEN}",
        "Content-Type": "application/json",
        "Accept": "application/vnd.github.v3+json",
        "X-GitHub-Api-Version": "2022-11-28",
    }

    response = requests.get(url, headers=headers)

    return response.json()


def get_commits_for_asset(owner, repo, path, since=None):
    base_url = "https://api.github.com"

    # Endpoint
    endpoint = f"/repos/{owner}/{repo}/commits"
    url = base_url + endpoint

    # Headers
    headers = {
        "Authorization": f"token {GITHUB_PERSONAL_ACCESS_TOKEN}",
        "Content-Type": "application/json",
        "Accept": "application/vnd.github.v3+json",
        "X-GitHub-Api-Version": "2022-11-28",
    }

    params = {
        "path": path,
        # validate date is in ISO 8601 format
        "since": since.isoformat() if isinstance(since, datetime) else since,
    }

    response = requests.get(url, headers=headers, params=params)

    return response.json()


def get_repo_contents(owner, repo, path, commit_sha=None):
    base_url = "https://api.github.com"

    # Endpoint
    endpoint = f"/repos/{owner}/{repo}/contents/{path}"
    url = base_url + endpoint

    # Headers
    headers = {
        "Authorization": f"token {GITHUB_PERSONAL_ACCESS_TOKEN}",
        "Content-Type": "application/json",
        "Accept": "application/vnd.github.v3+json",
        "X-GitHub-Api-Version": "2022-11-28",
    }

    params = {"ref": commit_sha}

    response = requests.get(url, headers=headers, params=params)

    # check status
    if response.status_code != 200:
        logger.debug(json.dumps(response.json()))
        raise Exception(f"Asset not found at path '{path}'")
    else:
        return response.json()


def get_repo_asset(owner, repo, path, commit_sha=None):
    """Download a file from a GitHub repository"""
    repo_contents = get_repo_contents(owner, repo, path, commit_sha)

    res = requests.get(repo_contents["download_url"])

    if res.status_code != 200:
        logger.error(f"Status code {res.status_code} for {path}")
        raise Exception(f"Error downloading {path}")

    return res.text


def get_branches(owner, repo):
    """Fetch branches for a GitHub repository"""

    base_url = "https://api.github.com"

    # Endpoint
    endpoint = f"/repos/{owner}/{repo}/branches"
    url = base_url + endpoint

    # Headers
    headers = {
        "Authorization": f"token {GITHUB_PERSONAL_ACCESS_TOKEN}",
        "Content-Type": "application/json",
        "Accept": "application/vnd.github.v3+json",
        "X-GitHub-Api-Version": "2022-11-28",
    }

    response = requests.get(url, headers=headers)
    return response.json()


def get_branch(owner, repo, branch_name):
    """Fetch branches for a GitHub repository"""

    base_url = "https://api.github.com"

    # Endpoint
    endpoint = f"/repos/{owner}/{repo}/branches/{branch_name}"
    url = base_url + endpoint

    # Headers
    headers = {
        "Authorization": f"token {GITHUB_PERSONAL_ACCESS_TOKEN}",
        "Content-Type": "application/json",
        "Accept": "application/vnd.github.v3+json",
        "X-GitHub-Api-Version": "2022-11-28",
    }

    response = requests.get(url, headers=headers)
    return response.json()


# Function to fetch pull requests
def get_pull_requests(owner, repo):
    url = f"https://api.github.com/repos/{owner}/{repo}/pulls?state=all"

    # Headers
    headers = {
        "Authorization": f"token {GITHUB_PERSONAL_ACCESS_TOKEN}",
        "Content-Type": "application/json",
        "Accept": "application/vnd.github.v3+json",
        "X-GitHub-Api-Version": "2022-11-28",
    }
    response = requests.get(url, headers=headers)

    if response.status_code == 200:
        return response.json()
    else:
        print(f"Error: {response.status_code}")
        return []


# def merge_release_version_with_commit(unique_shas, release_versions):
#     # Convert release_versions to a dictionary for easier lookup
#     release_versions_dict = {sha: version for sha, version in release_versions}

#     # Merge the arrays using a list comprehension
#     merged_data = [(release_versions_dict[sha], sha, date) for sha, date in unique_shas]

#     # sort by date
#     merged_data.sort(key=lambda x: x[2])

#     return merged_data


def flatten_json_records(records, skip_fields=[]):
    """Flatten a list of JSON records."""
    return [flatten_json(record, skip_fields=skip_fields) for record in records]


def select_keys(d, keys):
    """Selects keys from a dictionary"""
    return {k: v for k, v in d.items() if k in keys}


def select_fields(dataset, fields):
    """Select the fields for each record in an array of JSON objects"""
    return [select_keys(x, fields) for x in dataset]


def rename_keys(d, key_names_map):
    """Rename keys in a dictionary"""
    return {key_names_map[k]: v for k, v in d.items()}


def rename_fields(dataset: List[dict], key_names_map: dict[str, str]):
    """Rename fields in each record in an array of JSON objects"""
    return [rename_keys(x, key_names_map) for x in dataset]


def find_text(pattern, input_str):
    match = re.search(pattern, input_str)
    if match:
        text = match.group(0)
        return text


def get_release_version_for_commit(commit: dict, **kwargs) -> int:
    sha = commit["sha"]
    asset_path = kwargs["asset_path"]
    release_version_regex = kwargs["metadata_regex"]
    allele_list = get_repo_asset(
        GITHUB_REPOSITORY_OWNER, GITHUB_REPOSITORY_NAME, asset_path, sha
    )
    release_version = find_text(release_version_regex, allele_list)
    if release_version is None:
        raise Exception(f"Release version not found for commit {sha}")
    return int(release_version.replace(".", "")[:4])


def filter_nulls(items: List[Dict[str, str]]) -> List[Dict[str, str]]:
    return [x for x in items if x is not None]


def filter_null_fields(items: dict) -> dict:
    return { k: v for k, v in items.items() if v is not None }


def sort_execution_state_items(
    execution_state_items: List[Dict[str, str]], ascending=False
) -> List[Dict[str, str]]:
    return sorted(
        execution_state_items,
        key=lambda x: x["commit"].date_utc,
        reverse=(not ascending),
    )


def process_execution_state_item(
    commit: Dict[str, str],
    repository_config: RepositoryConfig,
    target_metadata_config: TargetMetadataConfig,
    limit: int = None,
) -> Dict[str, str]:
    errors = 0
    sha = commit["sha"]

    for config in target_metadata_config.items:
        try:
            logger.info(
                f"Getting release version for sha {sha} from {config.asset_path}"
            )
            release_version = get_release_version_for_commit(commit, **config.dict())
            logger.info(f"Found release version {release_version} ({sha})")

            result = {
                "repository": repository_config,
                "commit": Commit(**commit),
                "execution": ExecutionDetailsConfig(
                    version=release_version,
                    status="NOT_PROCESSED",
                    date_utc=None,
                    input_parameters=None,
                )
            }

        except Exception as e:
            # This is because Allelelist.txt for certain commits doesn't contain the release version or name
            # Need to find another file that indicates the release version should be small
            errors += 1
            logger.error(f"Error processing commit {sha}: {e}")
            if errors == len(target_metadata_config.items):
                # logger.error(f"Max errors reached. Exiting loop.")
                raise e
            else:
                continue

        # return error count and increment outside this function

        # TODO deserialize to ExecutionStateItem
        return result


def parallel_process_execution_state_items(
    commits: List[Dict[str, str]],
    repository_config: RepositoryConfig,
    target_metadata_config: TargetMetadataConfig,
    limit: int = None,
):
    execution_state_items = []
    num_cores = multiprocessing.cpu_count()
    num_threads = max(1, num_cores - 1)  # Reserve one core for other processes
    num_threads = min(6, num_cores)  # limit threads to avoid GitHub API rate limit

    # Create a ThreadPoolExecutor with the specified number of threads
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        # Submit the process_commit function for each commit to the executor
        futures = [
            executor.submit(
                process_execution_state_item, commit, repository_config, target_metadata_config
            )
            for commit in commits[:limit]
        ]

        # Collect the results as they complete
        execution_state_items = [future.result() for future in as_completed(futures)]

    return [
        ExecutionStateItem(**item)
        for item in sort_execution_state_items(filter_nulls(execution_state_items))
    ]


# limit is int or None
# @cache_pickle
def process_execution_state_items(
    commits: List[Dict[str, str]],
    repository_config: RepositoryConfig,
    target_metadata_config: TargetMetadataConfig,
    limit: None = None,
    parallel: str = False,
) -> List[Dict[str, str]]:
    if parallel == True:
        if limit:
            logger.warning("'limit' will not work if parallel processing is enabled")
        return parallel_process_execution_state_items(
            commits=commits,
            repository_config=repository_config,
            target_metadata_config=target_metadata_config,
            limit=limit,
        )
    else:
        execution_state_items = []
        for commit in commits[:limit]:
            execution_state_items.append(
                process_execution_state_item(
                    commit=commit,
                    repository_config=repository_config,
                    target_metadata_config=target_metadata_config,
                    limit=limit,
                )
            )

        return [
            ExecutionStateItem(**item)
            for item in sort_execution_state_items(filter_nulls(execution_state_items))
        ]
