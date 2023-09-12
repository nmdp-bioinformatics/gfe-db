"""
Builds the execution state for the given repository source from the static repository source configuration (`source-config.json`).
"""
import os
import sys

# for dev, local path to gfe-db modules
# ./gfe-db/pipeline/lambda_layers/gfe_db_models (use absolute path)
sys.path.append(os.environ["GFEDBMODELS_PATH"])

from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
# from datetime import datetime
import json
print(json.dumps(sys.path, indent=4))
# from gfedbmodels.constants import (
#     infra,
#     pipeline
# )
from gfedbmodels.utils import (
    get_utc_now,
    paginate_commits,
    select_fields,
    flatten_json_records,
    select_keys,
    rename_fields,
    filter_nested_nulls,
)
from gfedbmodels.types import (
    SourceConfig, 
    RepositoryConfig, 
    ExecutionState
)
from gfedbmodels.ingest import (
    process_execution_state_items
)

# Environment variables
GITHUB_REPOSITORY_OWNER = os.environ["GITHUB_REPOSITORY_OWNER"]
GITHUB_REPOSITORY_NAME = os.environ["GITHUB_REPOSITORY_NAME"]
GITHUB_PERSONAL_ACCESS_TOKEN = os.environ["GITHUB_PERSONAL_ACCESS_TOKEN"]
# DATA_BUCKET_NAME = infra.params.DataBucketName


if __name__ == "__main__":

    utc_now = get_utc_now()
    
    # Paths
    try:
        output_dir = Path(sys.argv[1])
    except IndexError:
        raise ValueError("Output directory must be specified as first argument")

    with open(output_dir / "source-config.json", "r") as f:
        source_config = SourceConfig(**json.load(f))

    # TODO FIX not returning all commits from repo, integrate pygethub
    # Fetch all commits from repo using GitHub API, will be cached
    logger.info("Fetching all commits from repo using GitHub API")
    all_commits = paginate_commits(GITHUB_REPOSITORY_OWNER, GITHUB_REPOSITORY_NAME, token=GITHUB_PERSONAL_ACCESS_TOKEN)

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
    target_metadata_config = source_config.repositories[
        GITHUB_REPOSITORY_OWNER + "/" + GITHUB_REPOSITORY_NAME
    ].target_metadata_config
    excluded_commit_shas = source_config.repositories[
        GITHUB_REPOSITORY_OWNER + "/" + GITHUB_REPOSITORY_NAME
    ].excluded_commit_shas.values
    commits = [
        commit for commit in commits if commit["sha"] not in excluded_commit_shas
    ]

    # Build ExecutionStateItem records from raw commits using thread pooling
    logger.info("Building execution state")
    execution_state_items = process_execution_state_items(
        timestamp=utc_now,
        commits=commits,
        # TODO remove default params from state table, they are retrieved from source config file in S3
        repository_config=RepositoryConfig(
            **select_keys(
                source_config.repositories[
                    GITHUB_REPOSITORY_OWNER + "/" + GITHUB_REPOSITORY_NAME
                ].model_dump(),
                ["owner", "name", "url"],
            )
        ),
        target_metadata_config=target_metadata_config,  # Infers release version from file contents
        token=GITHUB_PERSONAL_ACCESS_TOKEN,
        parallel=False,
    )

    # Package records as ExecutionState object to seed table
    execution_state = ExecutionState(
        **{
            "created_utc": utc_now,
            "items": execution_state_items,
        }
    )

    # Updates the source config file but does not actually build it
    source_config.created_utc, source_config.updated_utc = utc_now, utc_now

    logger.info(
        f"Writing execution state to {str(output_dir / 'execution-state.json')}"
    )

    # write ExecutionState locally
    with open(output_dir / "execution-state.json", "w") as f:
        json.dump(filter_nested_nulls(execution_state.model_dump()), f, indent=4)

    logger.info(f"Updating source config in {str(output_dir / 'source-config.json')}")

    # write SourceConfig locally
    with open(output_dir / f"source-config.json", "w") as f:
        json.dump(source_config.model_dump(), f, indent=4)
