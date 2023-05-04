"""
Builds the source config file for the given repository source.
"""
import os
import sys
sys.path.append(
    "/Users/ammon/Projects/nmdp-bioinformatics/02-Repositories/gfe-db/gfe-db/pipeline/functions/check_source_update"
)
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
from datetime import datetime

utc_now = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
import json

# these libraries are shared from the check_source_update function
from src.utils.types import (
    InputParameters,
    ExecutionStateItem,
    RepositoryConfig,
    SourceConfig,
    # SourceConfigBase
)
from src.utils import (
    read_source_config_base,
    paginate_commits,
    select_fields,
    flatten_json_records,
    select_keys,
    rename_fields,
    process_execution_state_items,
)

# Environment variables
APP_NAME = os.environ["APP_NAME"]
GITHUB_REPOSITORY_OWNER = "ANHIG"  # os.environ["GITHUB_REPOSITORY_OWNER"]
GITHUB_REPOSITORY_NAME = "IMGTHLA"  # os.environ["GITHUB_REPOSITORY_NAME"]
DATA_BUCKET_NAME = os.environ["DATA_BUCKET_NAME"]
PIPELINE_SOURCE_CONFIG_BASE_S3_PATH = os.environ["PIPELINE_SOURCE_CONFIG_BASE_S3_PATH"]


if __name__ == "__main__":

    # Get base source config
    source_config_base = read_source_config_base(DATA_BUCKET_NAME, PIPELINE_SOURCE_CONFIG_BASE_S3_PATH)

    # Paths
    # TODO arg
    output_dir = Path(f"{APP_NAME}/pipeline/config")

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
        "html_url": "html_url",
    }

    commits = rename_fields(all_commits_flat, commit_key_names_map)

    # Get the release version for each commit
    target_metadata_config = source_config_base.repositories[GITHUB_REPOSITORY_OWNER+"/"+GITHUB_REPOSITORY_NAME].target_metadata_config.values
    excluded_commit_shas = source_config_base.repositories[GITHUB_REPOSITORY_OWNER+"/"+GITHUB_REPOSITORY_NAME].excluded_commit_shas.values
    commits = [commit for commit in commits if commit["sha"] not in excluded_commit_shas]

    # Build ExecutionStateItem list using thread pooling    
    execution_state_items = process_execution_state_items(
        commits=commits,
        target_metadata_config=target_metadata_config,
        parallel=False,
    )

    # # convert utc_now string to datetime object using YYYYMMDDHHMM format
    # utc_now_version = datetime.strptime(utc_now, "%Y-%m-%dT%H:%M:%SZ").strftime(
    #     "%Y%m%d-%H%M"
    # )

    # # TODO get from argparse
    # # write SourceConfig locally
    # with open(output_dir / f"source-config.json", "w") as f:
    #     json.dump(source_config.dict(), f, indent=4)

exit(0)
