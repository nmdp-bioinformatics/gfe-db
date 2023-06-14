"""
Builds the execution state for the given repository source from the static repository source configuration (`source-config.json`).
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
    filter_nested_nulls
)

# Environment variables
APP_NAME = os.environ["APP_NAME"]
GITHUB_REPOSITORY_OWNER = os.environ["GITHUB_REPOSITORY_OWNER"]
GITHUB_REPOSITORY_NAME = os.environ["GITHUB_REPOSITORY_NAME"]
DATA_BUCKET_NAME = os.environ["DATA_BUCKET_NAME"]


if __name__ == "__main__":

    # Paths
    output_dir = Path(sys.argv[1])

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

    # rename fields for state table model
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

    # Build ExecutionStateItem records from raw commits using thread pooling    
    logger.info("Building execution state")
    execution_state_items = process_execution_state_items(
        timestamp=utc_now,
        commits=commits,
        # TODO remove default params from state table, they are retrieved from source config file in S3
        repository_config=RepositoryConfig(
            **select_keys(
                source_config.repositories[GITHUB_REPOSITORY_OWNER+"/"+GITHUB_REPOSITORY_NAME].dict(), 
                ["owner", "name", "url"]
            )
        ),
        target_metadata_config=target_metadata_config, # Infers release version from file contents
        parallel=True,
    )

    # Package records as ExecutionState object to seed table
    execution_state = ExecutionState(**{
        "created_utc": utc_now,
        "items": execution_state_items,
    })

    # Updates the source config file but does not actually build it
    source_config.created_utc, source_config.updated_utc = utc_now, utc_now

    logger.info(f"Writing execution state to {str(output_dir / 'execution-state.json')}")

    # write ExecutionState locally
    with open(output_dir / "execution-state.json", "w") as f:
        json.dump(filter_nested_nulls(execution_state.dict()), f, indent=4)

    logger.info(f"Updating source config in {str(output_dir / 'source-config.json')}")

    # write SourceConfig locally
    with open(output_dir / f"source-config.json", "w") as f:
        json.dump(source_config.dict(), f, indent=4)
