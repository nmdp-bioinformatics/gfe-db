"""
Builds the source config file for the given repository source.
"""
import os
import sys
sys.path.append("/Users/ammon/Projects/nmdp-bioinformatics/02-Repositories/gfe-db/scripts/")
from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
from typing import List, Dict
from datetime import datetime
utc_now = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
import re
import json
import multiprocessing
from concurrent.futures import ThreadPoolExecutor, as_completed
from build_source_config.utils.types import (
    Commit,
    InputParameters,
    ExecutionHistoryItem,
    RepositoryConfig,
    SourceConfig
)
from build_source_config.utils import (
    cache_pickle,
    paginate_commits,
    flatten_json,
    get_repo_asset
)

# Environment variables
APP_NAME = os.environ["APP_NAME"]
GITHUB_REPOSITORY_OWNER = "ANHIG" # os.environ["GITHUB_REPOSITORY_OWNER"]
GITHUB_REPOSITORY_NAME = "IMGTHLA" # os.environ["GITHUB_REPOSITORY_NAME"]
DATA_BUCKET_NAME = os.environ["DATA_BUCKET_NAME"]
PIPELINE_CONFIG_S3_PATH = os.environ["PIPELINE_CONFIG_S3_PATH"]

# Paths

# TODO arg
output_dir = Path(f"{APP_NAME}/pipeline/config")


def flatten_json_records(records):
    """Flatten a list of JSON records."""
    return [flatten_json(record) for record in records]


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
        sha = commit['sha']
        asset_name = kwargs['asset_name']
        release_version_regex = kwargs['release_version_regex']
        allele_list = get_repo_asset(GITHUB_REPOSITORY_OWNER, GITHUB_REPOSITORY_NAME, asset_name, sha)
        release_version = find_text(release_version_regex, allele_list)
        if release_version is None:
            raise Exception(f"Release version not found for commit {sha}")
        return int(release_version.replace(".", "")[:4])


def process_execution_history_item(commit: Dict[str, str], asset_configs: Dict[str, str], max_errors: int, errors: int) -> Dict[str, str]:

    for asset_config in asset_configs:
        try:
            logger.info(f"Trying config {json.dumps(asset_config)}")
            sha = commit['sha']
            date = commit['date_utc']
            logger.debug(f"Getting release version for sha {sha} and date {date}")

            release_version = get_release_version_for_commit(commit, **asset_config)
            logger.info(f"Found release version {release_version} for sha {sha}")

            result = {
                "version": release_version,
                "execution_date_utc": None,
                "commit": Commit(**commit),
                "input_parameters": None,
                "status": None
            }
        except Exception as e:
            # This is because Allelelist.txt for certain commits doesn't contain the release version or name
            # Need to find another file that indicates the release version should be small
            errors += 1
            logger.error(f"Error processing commit {commit['sha']}: {e}")
            if errors >= max_errors:
                logger.error(f"Max errors reached. Exiting loop.")
                break
            else:
                continue

        # return error count and increment outside this function
        return result

def filter_nulls(items: List[Dict[str, str]]) -> List[Dict[str, str]]:
    return [x for x in items if x is not None]


def sort_execution_history_items(execution_history_items: List[Dict[str, str]], ascending=False) -> List[Dict[str, str]]:
    return sorted(execution_history_items, key=lambda x: x.commit.date_utc, reverse=(not ascending))

# limit is int or None
@cache_pickle
def process_execution_history_items(commits: List[Dict[str, str]], asset_configs: List[Dict[str, str]], limit: int, max_errors: int) -> List[Dict[str, str]]:
    execution_history_items = []
    errors = 0
    num_cores = multiprocessing.cpu_count()
    num_threads = max(1, num_cores - 1)  # Reserve one core for other processes

    # Create a ThreadPoolExecutor with the specified number of threads
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        # Submit the process_commit function for each commit to the executor
        futures = [
            executor.submit(process_execution_history_item, commit, asset_configs, max_errors, errors)
            for commit in commits[:limit]
        ]

        # Collect the results as they complete
        execution_history_items = [future.result() for future in as_completed(futures)]

    return sort_execution_history_items(filter_nulls(execution_history_items))

if __name__ == "__main__":

    # Fetch all commits from repo using GitHub API, will be cached
    all_commits = paginate_commits(GITHUB_REPOSITORY_OWNER, GITHUB_REPOSITORY_NAME)
    
    # filter by chosen commit keys
    commit_keys = ["sha", "commit", "html_url"]
    all_commits = select_fields(all_commits, commit_keys)
    
    # flatten JSON records
    all_commits_flat = flatten_json_records(all_commits)
    commit_keys = ["sha", "commit.author.date", "commit.message", "html_url"]
    all_commits_flat = [select_keys(x, commit_keys) for x in all_commits_flat]

    # refactor to use rename_fields
    commit_key_names_map = {
        "sha": "sha", 
        "commit.author.date": "date_utc", 
        "commit.message": "message", 
        "html_url": "html_url"
    }

    commits = rename_fields(all_commits_flat, commit_key_names_map)

    # Get the release version for each commit    
    max_errors = 5
    limit = None
    asset_configs = [
        {
            "asset_name": "alignments/V_nuc.txt", # commits from 3a71348 to current
            "release_version_regex": r'[1-9]\d{0,1}\.[1-9]\d{0,2}\.0(?:\.\d{1,2})?(?=\s|$)'
        },
        {
            "asset_name": "aligments/V_nuc.txt", # commits from 8632b0d to 3645f26
            "release_version_regex": r'[1-9]\d{0,1}\.[1-9]\d{0,2}\.0(?:\.\d{1,2})?(?=\s|$)',
        },
        {
            "asset_name": "Alignments/V_nuc.txt", # commits from af54d28 to 9d8f585
            "release_version_regex": r'[1-9]\d{0,1}\.[1-9]\d{0,2}\.0(?:\.\d{1,2})?(?=\s|$)'
        },
        {
            "asset_name": "V_nuc.txt", # all commits before 08e0ef9
            "release_version_regex": r'[1-9]\d{0,1}\.[1-9]\d{0,2}\.0(?:\.\d{1,2})?(?=\s|$)',
        },
    ]

    # Build ExecutionHistoryItem list
    execution_history_items = process_execution_history_items(commits, asset_configs, limit, max_errors)

    # Build RepositoryConfig
    repository_path = f"{GITHUB_REPOSITORY_OWNER}/{GITHUB_REPOSITORY_NAME}"
    tracked_assets = ["hla.dat", "msf/"]

    source_config = SourceConfig(**{
        "created_at_utc": utc_now,
        "updated_at_utc": utc_now,
        "repositories": {
            repository_path: RepositoryConfig(**{
                "owner": GITHUB_REPOSITORY_OWNER,
                "name": GITHUB_REPOSITORY_NAME,
                "url": f"https://github.com/{repository_path}",
                "tracked_assets": tracked_assets,
                # TODO fetch default input parameters from S3 so they are easy to update 
                "default_input_parameters": InputParameters(**{
                    "align": "False",
                    "kir": "False",
                    "mem_profile": "False",
                    "limit": "1000",
                }),
                "execution_history": [ ExecutionHistoryItem(**item) for item in execution_history_items ]
            })
        }
    })

    # convert utc_now string to datetime object using YYYYMMDDHHMM format
    utc_now_version = datetime.strptime(utc_now, "%Y-%m-%dT%H:%M:%SZ").strftime("%Y%m%d-%H%M")

    # write SourceConfig locally
    with open(output_dir / f"source-config-v{utc_now_version}.json", "w") as f:
        json.dump(source_config.dict(), f, indent=4)

exit(0)