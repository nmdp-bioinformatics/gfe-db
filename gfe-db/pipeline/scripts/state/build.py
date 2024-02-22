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
from pygethub import list_commits, list_branches, GitHubPaginator
from gfedbmodels.utils import (
    get_utc_now,
    select_fields,
    flatten_json_records,
    select_keys,
    rename_fields,
    filter_nested_nulls,
    get_commit,
    cache_pickle
)
from gfedbmodels.types import (
    SourceConfig, 
    RepositoryConfig,
    Commit,
    ExecutionStateItem,
    ExecutionDetailsConfig,
    ExecutionState
)
from gfedbmodels.ingest import (
    process_execution_state_items,
    get_commits_by_branch,
    sort_execution_state_items
)

# Environment variables
GITHUB_REPOSITORY_OWNER = os.environ["GITHUB_REPOSITORY_OWNER"]
GITHUB_REPOSITORY_NAME = os.environ["GITHUB_REPOSITORY_NAME"]
GITHUB_PERSONAL_ACCESS_TOKEN = os.environ["GITHUB_PERSONAL_ACCESS_TOKEN"]
# DATA_BUCKET_NAME = infra.params.DataBucketName

# @cache_pickle
def get_branch_commits(branches):

    # For each entry in all-branches, get the commit data and build the execution state item
    execution_state_items = []
    for item in branches:

        if item['name'].lower() == 'latest':
            continue
        
        release_version = item['name']
        sha = item['commit']['sha']

        commit_json = get_commit(GITHUB_REPOSITORY_OWNER, GITHUB_REPOSITORY_NAME, GITHUB_PERSONAL_ACCESS_TOKEN, sha)
        assert sha == commit_json['sha']

        logger.info(f'Retrieving data for {sha}')
        
        commit = Commit(
            sha=commit_json['sha'],
            date_utc=commit_json['commit']['author']['date'],
            message=commit_json['commit']['message'],
            html_url=commit_json['html_url']
        )

        execution_state_item = ExecutionStateItem(
            created_utc=utc_now,
            updated_utc=utc_now,
            commit=commit,
            execution=ExecutionDetailsConfig(
                version=release_version,
                status="NOT_PROCESSED",
                date_utc=None,
                input_parameters=None,
            ),
            # error=None,
            # s3_path=None,
            repository=RepositoryConfig(
                **select_keys(
                    source_config.repositories[
                        GITHUB_REPOSITORY_OWNER + "/" + GITHUB_REPOSITORY_NAME
                    ].model_dump(),
                    ["owner", "name", "url"],
                )
            ),
        )
        execution_state_items.append(execution_state_item)

    return execution_state_items


def build_execution_state(branches, utc_now=None):

    utc_now = utc_now or get_utc_now()

    # Create ExecutionStateItems array from branch/commit sha pairs
    execution_state_items = get_branch_commits(branches)

    # Sort execution state items by date descending
    execution_state_items = sorted(
        execution_state_items,
        key=lambda x: x.commit.date_utc,
        reverse=True
    )

    # Package records as ExecutionState object to seed table
    execution_state = ExecutionState(
        **{
            "created_utc": utc_now,
            "items": execution_state_items,
        }
    )

    return execution_state

if __name__ == "__main__":

    utc_now = get_utc_now()
    
    # Paths
    try:
        output_dir = Path(sys.argv[1])
    except IndexError:
        raise ValueError("Output directory must be specified as first argument")

    with open(output_dir / "source-config.json", "r") as f:
        source_config = SourceConfig(**json.load(f))

    # Fetch all commits from repo using GitHub API, will be cached
    logger.info("Processing source repository data")

    # TODO add requests session for user-agent tracking
    paginator = GitHubPaginator(GITHUB_PERSONAL_ACCESS_TOKEN)

    ### COMMITS BY BRANCHES ###
    branch_pages = paginator.get_paginator(
        list_branches, 
        owner=GITHUB_REPOSITORY_OWNER, 
        repo=GITHUB_REPOSITORY_NAME,
        user_agent="nmdp-bioinformatics-gfe-db-state-builder/1.0"
    )
    all_branches = list(branch_pages)

    # # extract the branch names
    # branch_names = [branch["name"] for branch in all_branches]

    # # For each entry in all-branches, get the commit data and build the execution state item
    # execution_state_items = []
    # for item in all_branches:

    #     if item['name'].lower() == 'latest':
    #         continue
        
    #     release_version = item['name']
    #     sha = item['commit']['sha']

    #     commit_json = get_commit(GITHUB_REPOSITORY_OWNER, GITHUB_REPOSITORY_NAME, GITHUB_PERSONAL_ACCESS_TOKEN, sha)
    #     assert sha == commit_json['sha']

    #     logger.info(f'Retrieving data for {sha}')
        
    #     commit = Commit(
    #         sha=commit_json['sha'],
    #         date_utc=commit_json['commit']['author']['date'],
    #         message=commit_json['commit']['message'],
    #         html_url=commit_json['html_url']
    #     )

    #     execution_state_item = ExecutionStateItem(
    #         created_utc=utc_now,
    #         updated_utc=utc_now,
    #         commit=commit,
    #         execution=ExecutionDetailsConfig(
    #                     version=release_version,
    #                     status="NOT_PROCESSED",
    #                     date_utc=None,
    #                     input_parameters=None,
    #         ),
    #         # error=None,
    #         # s3_path=None,
    #         repository=RepositoryConfig(
    #             **select_keys(
    #                 source_config.repositories[
    #                     GITHUB_REPOSITORY_OWNER + "/" + GITHUB_REPOSITORY_NAME
    #                 ].model_dump(),
    #                 ["owner", "name", "url"],
    #             )
    #         ),
    #     )
    #     execution_state_items.append(execution_state_item)

    # Create ExecutionStateItems object from branch/commit sha pairs
    # Commit fields
    """
    class Commit(BaseModel):
        sha: str
        date_utc: str
        message: Optional[str] = None
        html_url: str
    """
    # ExecutionStateItems fields
    """
    class ExecutionStateItem(BaseModel):
        created_utc: Optional[str] = None
        updated_utc: Optional[str] = None
        repository: Optional[RepositoryConfig]
        commit: Commit
        execution: ExecutionDetailsConfig
        error: Optional[ExecutionError] = None
        s3_path: Optional[str] = None
    """

    # # COMMITS
    # commit_pages = paginator.get_paginator(
    #     list_commits, 
    #     owner=GITHUB_REPOSITORY_OWNER, 
    #     repo=GITHUB_REPOSITORY_NAME, 
    #     user_agent="nmdp-bioinformatics-gfe-db-state-builder/1.0")
    # all_commits = list(commit_pages)

    # # filter by chosen commit keys
    # commit_keys = ["sha", "commit", "html_url"]
    # all_commits = select_fields(all_commits, commit_keys)

    # # flatten JSON records
    # all_commits_flat = flatten_json_records(all_commits)
    # commit_keys = ["sha", "commit.author.date", "commit.message", "html_url"]
    # all_commits_flat = [select_keys(x, commit_keys) for x in all_commits_flat]

    # # rename fields for state table model
    # commit_key_names_map = {
    #     "sha": "sha",
    #     "commit.author.date": "date_utc",
    #     "commit.message": "message",
    #     "html_url": "html_url",
    # }

    # commits = rename_fields(all_commits_flat, commit_key_names_map)

    # # Get the release version for each commit
    # target_metadata_config = source_config.repositories[
    #     GITHUB_REPOSITORY_OWNER + "/" + GITHUB_REPOSITORY_NAME
    # ].target_metadata_config

    # excluded_commit_shas = source_config.repositories[
    #     GITHUB_REPOSITORY_OWNER + "/" + GITHUB_REPOSITORY_NAME
    # ].excluded_commit_shas.values

    # commits = [
    #     commit for commit in commits if commit["sha"] not in excluded_commit_shas
    # ]

    # # Build ExecutionStateItem records from raw commits using thread pooling
    # logger.info("Building execution state")
    # error_shas, execution_state_items = process_execution_state_items(
    #     timestamp=utc_now,
    #     commits=commits,
    #     # TODO remove default params from state table, they are retrieved from source config file in S3
    #     repository_config=RepositoryConfig(
    #         **select_keys(
    #             source_config.repositories[
    #                 GITHUB_REPOSITORY_OWNER + "/" + GITHUB_REPOSITORY_NAME
    #             ].model_dump(),
    #             ["owner", "name", "url"],
    #         )
    #     ),
    #     target_metadata_config=target_metadata_config,  # Infers release version from file contents
    #     token=GITHUB_PERSONAL_ACCESS_TOKEN,
    #     parallel=True,
    # )

    # # Sort execution state items by date descending
    # execution_state_items = sorted(
    #     execution_state_items,
    #     key=lambda x: x.commit.date_utc,
    #     reverse=True
    # )

    # # Package records as ExecutionState object to seed table
    # execution_state = ExecutionState(
    #     **{
    #         "created_utc": utc_now,
    #         "items": execution_state_items,
    #     }
    # )

    execution_state = build_execution_state(all_branches, utc_now)

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

    logger.info("Execution state and source config updated")

    # # write error shas
    # with open(output_dir / "error-shas.json", "w") as f:
    #     json.dump(error_shas, f, indent=4)
