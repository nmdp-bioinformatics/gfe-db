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
from pygethub import list_branches, GitHubPaginator
from gfedbmodels.utils import get_utc_now, select_keys, filter_nested_nulls, get_commit
from gfedbmodels.types import (
    SourceConfig,
    RepositoryConfig,
    Commit,
    ExecutionStateItem,
    ExecutionDetailsConfig,
    ExecutionState,
)

# Environment variables
GITHUB_REPOSITORY_OWNER = os.environ["GITHUB_REPOSITORY_OWNER"]
GITHUB_REPOSITORY_NAME = os.environ["GITHUB_REPOSITORY_NAME"]
GITHUB_PERSONAL_ACCESS_TOKEN = os.environ["GITHUB_PERSONAL_ACCESS_TOKEN"]


def get_branch_commits(branches):

    # For each entry in all-branches, get the commit data and build the execution state item
    execution_state_items = []

    for item in branches:

        if item["name"].lower() == "latest":
            continue

        release_version = item["name"]
        sha = item["commit"]["sha"]

        logger.info(f"Retrieving data for {sha}")
        commit_json = get_commit(
            GITHUB_REPOSITORY_OWNER,
            GITHUB_REPOSITORY_NAME,
            GITHUB_PERSONAL_ACCESS_TOKEN,
            sha,
        )
        assert sha == commit_json["sha"]

        commit = Commit(
            sha=commit_json["sha"],
            date_utc=commit_json["commit"]["author"]["date"],
            message=commit_json["commit"]["message"],
            html_url=commit_json["html_url"],
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
        execution_state_items, key=lambda x: x.commit.date_utc, reverse=True
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

    paginator = GitHubPaginator(GITHUB_PERSONAL_ACCESS_TOKEN)

    ### COMMITS BY BRANCHES ###
    branch_pages = paginator.get_paginator(
        list_branches,
        owner=GITHUB_REPOSITORY_OWNER,
        repo=GITHUB_REPOSITORY_NAME,
        user_agent="nmdp-bioinformatics-gfe-db-state-builder/1.0",
    )
    all_branches = list(branch_pages)

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
