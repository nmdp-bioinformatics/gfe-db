import os
import logging
from typing import List, Dict, Union
from pydantic import BaseModel
from pathlib import Path
from itertools import chain, starmap
from datetime import datetime
import json
import pickle
import re
import requests
import multiprocessing
from concurrent.futures import ThreadPoolExecutor, as_completed
from botocore.exceptions import ClientError
from .types import (
    SourceConfig,
    RepositoryConfig,
    TargetMetadataConfig,
    Commit,
    ExecutionStateItem,
    ExecutionDetailsConfig,
)
from .constants import (
    infra,
    pipeline,
)

# Logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

AWS_REGION = os.environ["AWS_REGION"]

# TODO can call these directly in the functions instead of decarling and passing them in, they should be cached
GITHUB_PERSONAL_ACCESS_TOKEN = infra.secrets.GitHubPersonalAccessToken
GITHUB_REPOSITORY_OWNER = pipeline.params.GitHubSourceRepository["owner"]
GITHUB_REPOSITORY_NAME = pipeline.params.GitHubSourceRepository["name"]

cache_dir = Path(__file__).parent / "_cache"

# TODO clear cache
# TODO disable/enable cache for testing

def save_json_to_cache(data, var_name):
    """Saves data to cache directory"""
    if not cache_dir.exists():
        cache_dir.mkdir()

    # Handles different types of JSON representations, if it fails to serialize, remove the file
    try:
        if isinstance(data, dict) or isinstance(data, List):
            try:
                with open(cache_dir / var_name, "w") as f:
                    json.dump(data, f, indent=4)
            except:
                # assume it's a list of pydantic models
                with open(cache_dir / var_name, "w") as f:
                    json.dump([item.dict() for item in data], f, indent=4)
    except Exception as e:
        logger.error(f"Failed to serialize {var_name} to JSON: {e}")
        # remove the file if it exists
        if (cache_dir / var_name).exists():
            (cache_dir / var_name).unlink()


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


def flatten_json(data, sep=".", skip_fields=[], select_fields=[]):
    """Flatten a nested json file. For a list of dictionaries, use this
    inside a for loop before converting to pandas DataFrame.

    Args:
        data (dict): nested json file
        sep (str, optional): separator for flattened keys. Defaults to ".".
        skip_fields (list, optional): list of fields to skip. Defaults to [].
        select_fields (list, optional): list of output fields to select including the separator. Defaults to [].
    """

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
        # Keep unpacking the json file until all values are atomic elements (not data or list)
        data = dict(chain.from_iterable(starmap(unpack, data.items())))
        # Terminate condition: not any value in the json file is data or list
        if not any(isinstance(value, dict) for value in data.values()) and not any(
            isinstance(value, list) for value in data.values()
        ):
            break

    if len(skip_fields) > 0:
        data = {k: v for k, v in data.items() if k not in skip_fields}

    if len(select_fields) > 0:
        data = {k: v for k, v in data.items() if k in select_fields}

    return data


def read_s3_json(s3_client, bucket, key):
    """Reads config file containing the current state of branches in
    a GitHub repo"""

    try:
        response = s3_client.get_object(Bucket=bucket, Key=key)
        return json.loads(response["Body"].read().decode())

    except ClientError as err:
        logger.error(f"Failed to read config file to s3://{bucket}/{key}")
        raise err


def write_s3_json(s3_client, bucket, key, data):
    """Writes config file containing the current state of branches in
    a GitHub repo"""

    try:
        response = s3_client.put_object(Bucket=bucket, Key=key, Body=json.dumps(data).encode())

    except Exception as err:
        logger.error(
            f'Failed to write config file to s3://{bucket}/{key}. HTTPStatusCode: {response["ResponseMetadata"]["HTTPStatusCode"]}'
        )
        raise err


def read_source_config(s3_client, bucket, key):
    data = read_s3_json(s3_client, bucket, key)
    return SourceConfig(**data)


# def write_source_config(bucket, key, source_config: SourceConfig):
#     write_s3_json(bucket, key, source_config.dict())


def list_commits(owner, repo, **params):
    """Return a list of GitHub commits for the specified repository"""

    base_url = "https://api.github.com"

    # Endpoint
    endpoint = f"/repos/{owner}/{repo}/commits"

    url = base_url + endpoint

    # params = {
    #     "per_page": kwargs.get("per_page"),
    #     "page": kwargs.get("page"),
    # }

    # Headers
    headers = {
        "Authorization": f"token {GITHUB_PERSONAL_ACCESS_TOKEN}",
        "Content-Type": "application/json",
        "Accept": "application/vnd.github.v3+json",
        "X-GitHub-Api-Version": "2022-11-28",
    }

    response = requests.get(url, headers=headers, params=params)
    response.raise_for_status()

    return response.json()


@cache_json
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
    response.raise_for_status()

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
    response.raise_for_status()

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
    response.raise_for_status()

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
    response.raise_for_status()

    # # check status
    # if response.status_code != 200:
    #     logger.debug(json.dumps(response.json()))
    #     raise Exception(f"Asset not found at path '{path}'")
    # else:
    return response.json()


def get_repo_asset(owner, repo, path, commit_sha=None):
    """Download a file from a GitHub repository"""
    repo_contents = get_repo_contents(owner, repo, path, commit_sha)

    response = requests.get(repo_contents["download_url"])
    response.raise_for_status()

    # if response.status_code != 200:
    #     logger.error(f"Status code {response.status_code} for {path}")
    #     raise Exception(f"Error downloading {path}")

    return response.text


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
    response.raise_for_status()

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
    response.raise_for_status()

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


def flatten_json_records(
    data, sep=".", skip_fields=[], select_fields=[], filter_nulls=True
):
    """Flatten a list of JSON records."""
    if filter_nulls:
        return [
            filter_null_fields(
                flatten_json(
                    data=record,
                    sep=sep,
                    skip_fields=skip_fields,
                    select_fields=select_fields,
                )
            )
            for record in data
        ]
    else:
        return [
            flatten_json(
                data=record,
                sep=sep,
                skip_fields=skip_fields,
                select_fields=select_fields,
            )
            for record in data
        ]


def restore_nested_json(data: dict, split_on="."):
    """Restores a previously flattened JSON object into a nested JSON object.

    Args:
        data (dict): A flattened JSON object.

    Returns:
        dict: A nested JSON object.
    """
    result = {}
    for key, value in data.items():
        parts = key.split(split_on)
        current = result
        for part in parts[:-1]:
            if part not in current:
                current[part] = {}
            current = current[part]
        current[parts[-1]] = value
    return result


def find_text(pattern, input_str):
    match = re.search(pattern, input_str)
    if match:
        text = match.group(0)
        return text


def get_release_version_for_commit(commit: Union[Commit, dict], **kwargs) -> int:
    try:
        sha = commit["sha"]
    except:
        sha = commit.sha
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
    """Filter out null items from a list of dictionaries

    Args:
        items (List[Dict[str, str]]): A list of dictionaries

    Returns:
        List[Dict[str, str]]: A list of dictionaries with null items removed
    """
    return [x for x in items if x is not None]


def filter_null_fields(items: dict) -> dict:
    """Filter out null fields from a dictionary

    Args:
        items (dict): A dictionary

    Returns:
        dict: A dictionary with null fields removed
    """
    return {k: v for k, v in items.items() if v is not None}


def filter_nested_nulls(data: Union[dict, list]):
    """Filter out null fields from a nested dictionary or list of dictionaries

    Args:
        data (Union[dict, list]): A nested dictionary or list of dictionaries

    Returns:
        Union[dict, list]: A nested dictionary or list of dictionaries with null fields removed
    """
    if isinstance(data, list):
        return filter_nulls([filter_nested_nulls(i) for i in data])
    elif isinstance(data, dict):
        return filter_null_fields({k: filter_nested_nulls(v) for k, v in data.items()})
    else:
        return data


def sort_execution_state_items(
    execution_state_items: List[Dict[str, str]], ascending=False
) -> List[Dict[str, str]]:
    return sorted(
        execution_state_items,
        key=lambda x: x["commit"].date_utc,
        reverse=(not ascending),
    )


def process_execution_state_item(
    timestamp: str,
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
                "created_utc": timestamp,
                "repository": repository_config,
                "commit": Commit(**commit),
                "execution": ExecutionDetailsConfig(
                    version=release_version,
                    status="NOT_PROCESSED",
                    date_utc=None,
                    input_parameters=None,
                ),
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
    timestamp: str,
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
                process_execution_state_item,
                timestamp,
                commit,
                repository_config,
                target_metadata_config,
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
@cache_pickle
def process_execution_state_items(
    timestamp: str,
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
            timestamp=timestamp,
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
