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
from datetime import datetime
utc_now = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
import re
import json
from build_source_config.utils.types import SourceConfig
from build_source_config.utils import (
    cache,
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


def rename_fields(dataset, key_names_map):
    """Rename fields in each record in an array of JSON objects"""
    return [rename_keys(x, key_names_map) for x in dataset]


if __name__ == "__main__":
    # Fetch all commits from repo using GitHub API, will be cached
    all_commits = paginate_commits(GITHUB_REPOSITORY_OWNER, GITHUB_REPOSITORY_NAME)
    
    # filter by select_keys
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
    # TODO use multithreading to speed up
    # next we get the release version for each commit
    release_version_re = r"# version: IPD-IMGT/HLA (\d+\.\d+\.\d+)"
    execution_history_items = []
    errors = 0
    max_errors = 5
    limit = 5
    for idx, commit in enumerate(commits):
        try:
            sha = commit['sha']
            date = commit['date_utc']
            logger.debug(f"Getting release version for sha {sha} and date {date}")
            allele_list = get_repo_asset(GITHUB_REPOSITORY_OWNER, GITHUB_REPOSITORY_NAME, "Allelelist.txt", sha)
            release_version = int(re.search(release_version_re, allele_list).group(1).replace(".", ""))
            logger.info(f"Found release version {release_version} for sha {sha}")

            execution_history_items.append({
                "version": release_version,
                "execution_date_utc": None,
                "commit": commit,
                "input_parameters": None,
                "status": None
            })
        except Exception as e:
            errors += 1
            logger.error(f"Error processing commit {commit['sha']}: {e}")
            if errors >= max_errors:
                logger.error(f"Max errors reached. Exiting loop.")
                break

        if limit is not None:
            if idx+1 == limit:
                break

    # Build RepositoryConfig
    repository_path = f"{GITHUB_REPOSITORY_OWNER}/{GITHUB_REPOSITORY_NAME}"
    tracked_assets = ["hla.dat", "msf/"]

    base_source_config = {
        "created_at_utc": utc_now,
        "updated_at_utc": utc_now,
        "repositories": {
            repository_path: {
                "owner": GITHUB_REPOSITORY_OWNER,
                "name": GITHUB_REPOSITORY_NAME,
                "url": f"https://github.com/{repository_path}",
                "tracked_assets": tracked_assets,
                "default_input_parameters": {
                    "align": "False",
                    "kir": "False",
                    "mem_profile": "False",
                    "limit": "1000"
                }
            }
        }
    }

    # Deserialize SourceConfig
    base_source_config["repositories"][repository_path]["execution_history"] = execution_history_items
    source_config = SourceConfig(**base_source_config)

    # convert utc_now string to datetime object using YYYYMMDDHHMM format
    utc_now_version = datetime.strptime(utc_now, "%Y-%m-%dT%H:%M:%SZ").strftime("%Y%m%d-%H%M")


    # write SourceConfig locally
    with open(output_dir / f"source-config-v{utc_now_version}.json", "w") as f:
        json.dump(source_config.dict(), f, indent=4)

exit(0)