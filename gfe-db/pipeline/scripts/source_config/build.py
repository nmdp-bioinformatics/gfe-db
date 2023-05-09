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
    SourceConfig,
    RepositoryConfig,
    ExecutionState
)
from src.utils import (
    paginate_commits,
    select_fields,
    flatten_json_records,
    select_keys,
    rename_fields,
    process_execution_state_items,
)

# Environment variables
APP_NAME = os.environ["APP_NAME"]
GITHUB_REPOSITORY_OWNER = os.environ["GITHUB_REPOSITORY_OWNER"]
GITHUB_REPOSITORY_NAME = os.environ["GITHUB_REPOSITORY_NAME"]
DATA_BUCKET_NAME = os.environ["DATA_BUCKET_NAME"]
# PIPELINE_SOURCE_CONFIG_S3_PATH = os.environ["PIPELINE_SOURCE_CONFIG_S3_PATH"]


if __name__ == "__main__":

    # Paths
    # TODO arg
    output_dir = Path(f"{APP_NAME}/pipeline/config")

    # Get base source config
    # try:
    #     source_config = SourceConfig(**read_s3_json(DATA_BUCKET_NAME, PIPELINE_SOURCE_CONFIG_S3_PATH))
    # except Exception as e:
    #     logger.error(f"Error loading source config from {PIPELINE_SOURCE_CONFIG_S3_PATH}")
    #     logger.info("Trying local source config")
    with open(output_dir / "source-config.json", "r") as f:
        source_config = SourceConfig(**json.load(f))

    # Fetch all commits from repo using GitHub API, will be cached
    logger.info("Fetching all commits from repo using GitHub API")
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
    target_metadata_config = source_config.repositories[GITHUB_REPOSITORY_OWNER+"/"+GITHUB_REPOSITORY_NAME].target_metadata_config
    excluded_commit_shas = source_config.repositories[GITHUB_REPOSITORY_OWNER+"/"+GITHUB_REPOSITORY_NAME].excluded_commit_shas.values
    commits = [commit for commit in commits if commit["sha"] not in excluded_commit_shas]

    # Build ExecutionStateItem list using thread pooling    
    logger.info("Building execution state")
    execution_state_items = process_execution_state_items(
        commits=commits,
        repository_config=RepositoryConfig(
            **select_keys(
                source_config.repositories[GITHUB_REPOSITORY_OWNER+"/"+GITHUB_REPOSITORY_NAME].dict(), 
                ["owner", "name", "url", "default_input_parameters"]
            )
        ),
        target_metadata_config=target_metadata_config,
        parallel=True,
    )

    execution_state = ExecutionState(**{
        "created_utc": utc_now,
        "updated_utc": utc_now,
        "items": execution_state_items,
    })

    source_config.created_utc, source_config.updated_utc = utc_now, utc_now

    logger.info("Writing execution state to local file system")
    # TODO get output_dir from argparse
    # write ExecutionState locally
    with open(output_dir / f"execution-state.json", "w") as f:
        json.dump(execution_state.dict(), f, indent=4)

    logger.info("Writing source config to local file system")
    # TODO get output_dir from argparse
    # write SourceConfig locally
    with open(output_dir / f"source-config.json", "w") as f:
        json.dump(source_config.dict(), f, indent=4)
