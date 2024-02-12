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
    paginate_commits,
    select_fields,
    flatten_json_records,
    select_keys,
    rename_fields,
    filter_nested_nulls,
    cache_json
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

    # Fetch all commits from repo using GitHub API, will be cached
    logger.info("Processing source repository data")

    # TODO add requests session for user-agent tracking
    paginator = GitHubPaginator(GITHUB_PERSONAL_ACCESS_TOKEN)

    ### COMMITS BY BRANCHES ###
    branch_pages = paginator.get_paginator(list_branches, owner=GITHUB_REPOSITORY_OWNER, repo=GITHUB_REPOSITORY_NAME)
    all_branches = list(branch_pages)

    # extract the branch names
    branch_names = [branch["name"] for branch in all_branches]

    @cache_json
    def get_commits_by_branch(branch_name):
        
        # Fetch the commits for each branch
        commits_by_branch = {}
        for branch in branch_names:
            logger.info(f"Retrieving commits for branch {branch}")
            list_commits_params = {
                "owner": GITHUB_REPOSITORY_OWNER,
                "repo": GITHUB_REPOSITORY_NAME,
                "sha": branch,
            }
            branch_commit_pages = paginator.get_paginator(
                list_commits, 
                **list_commits_params,
                user_agent="nmdp-bioinformatics-gfe-db-state-builder/1.0")
            commits_by_branch[branch] = list(branch_commit_pages)

        # Create an array of all unique commits (using set method) in commits_by_branch and omit the branch information
        all_commits = set()
        for release, commits in commits_by_branch.items():
            all_commits.update([json.dumps(commit) for commit in commits])

        # covert back to dict
        all_commits = [json.loads(commit) for commit in all_commits]

        return all_commits

    # get all commits by branch
    all_commits = get_commits_by_branch(branch_names)

    # # COMMITS
    # commit_pages = paginator.get_paginator(
    #     list_commits, 
    #     owner=GITHUB_REPOSITORY_OWNER, 
    #     repo=GITHUB_REPOSITORY_NAME, 
    #     user_agent="nmdp-bioinformatics-gfe-db-state-builder/1.0")
    # all_commits = list(commit_pages)


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
    error_shas, execution_state_items = process_execution_state_items(
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
        parallel=True,
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

    logger.info("Execution state and source config updated")

    # write error shas
    with open(output_dir / "error-shas.json", "w") as f:
        json.dump(error_shas, f, indent=4)
