"""
Checks a GitHub repository for new commits and triggers data ingestion
"""
import os
import logging
import json
import boto3
from utils.constants import (
    GITHUB_REPOSITORY_OWNER,
    GITHUB_REPOSITORY_NAME,
)
from utils.types import ExecutionStateItem, Commit, to_datetime
from utils.utils import list_commits, process_execution_state_items

# set up logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Environment
AWS_REGION = os.environ["AWS_REGION"]
APP_NAME = os.environ["APP_NAME"]
STAGE = os.environ["STAGE"]

# AWS clients
ssm = boto3.client("ssm")
dynamodb = boto3.resource("dynamodb", region_name=AWS_REGION)

# Get SSM Parameters
table_name = ssm.get_parameter(
    Name=f'/{os.environ["APP_NAME"]}/{os.environ["STAGE"]}/{os.environ["AWS_REGION"]}/ExecutionStateTableName'
)["Parameter"]["Value"]

state_machine_arn = ssm.get_parameter(
    Name=f'/{os.environ["APP_NAME"]}/{os.environ["STAGE"]}/{os.environ["AWS_REGION"]}/UpdatePipelineArn'
)["Parameter"]["Value"]


def lambda_handler(event, context):
    logger.info(json.dumps(event))

    ### Retrieve Execution State ###
    # instantiate table
    table = dynamodb.Table(table_name)

    # get all the items from the table
    items = table.scan()["Items"]

    # sort items by commit.date_utc
    items = sorted(items, key=lambda x: x["commit.date_utc"], reverse=True)

    # Deserialize the items
    execution_state = [
        ExecutionStateItem(
            **{
                "version": item["version"],
                "commit": {
                    "sha": item["commit.sha"],
                    "message": item["commit.message"],
                    "date_utc": item["commit.date_utc"],
                    "html_url": item["commit.html_url"],
                },
                "status": item["status"],
                "execution_date_utc": item["execution_date_utc"],
                "input_parameters": item["input_parameters"],
            }
        )
        for item in items
    ]
    del items

    ### Retrieve Repository Commits ###
    # get the most recent 100 commits from github
    commits = list_commits(
        GITHUB_REPOSITORY_OWNER, GITHUB_REPOSITORY_NAME, per_page=100
    )
    repository_state = [
        Commit(
            **{
                "sha": commit["sha"],
                "message": commit["commit"]["message"],
                "date_utc": commit["commit"]["author"]["date"],
                "html_url": commit["html_url"],
            }
        )
        for commit in commits
    ]

    # get the newest commits where date_utc is more recent
    new_commits = [
        commit
        for commit in repository_state
        if commit.date_utc > execution_state[0].commit.date_utc
    ]

    ### Process New Commits ### 
    # TODO Note: current logic is making the assumption that commits retrieved from state have already been processed and are marked accordingly
    # TODO All commits from state should be be processed if status is null, because this indicates they have not been handled yet
    # if there are new commits, get the release version for each commit and deserialize
    if len(new_commits) > 0:
        asset_configs = [
            {
                "asset_name": "alignments/V_nuc.txt",  # commits from 3a71348 to current
                "release_version_regex": r"[1-9]\d{0,1}\.[1-9]\d{0,2}\.0(?:\.\d{1,2})?(?=\s|$)",
            },
            {
                "asset_name": "aligments/V_nuc.txt",  # commits from 8632b0d to 3645f26
                "release_version_regex": r"[1-9]\d{0,1}\.[1-9]\d{0,2}\.0(?:\.\d{1,2})?(?=\s|$)",
            },
            {
                "asset_name": "Alignments/V_nuc.txt",  # commits from af54d28 to 9d8f585
                "release_version_regex": r"[1-9]\d{0,1}\.[1-9]\d{0,2}\.0(?:\.\d{1,2})?(?=\s|$)",
            },
            {
                "asset_name": "V_nuc.txt",  # all commits before 08e0ef9
                "release_version_regex": r"[1-9]\d{0,1}\.[1-9]\d{0,2}\.0(?:\.\d{1,2})?(?=\s|$)",
            },
        ]

        # get the release version for each commit
        new_execution_state_items = process_execution_state_items(
            commits=[commit.dict() for commit in new_commits],
            asset_configs=asset_configs,
            limit=None,
            parallel=False,
        )

        # TODO move inside new_execution_state_items
        new_execution_state_items = [
            ExecutionStateItem(**item) for item in new_execution_state_items
        ]

        # group by release version and get most recent by commit date (max date_utc)
        unique_new_releases = list(set([item.version for item in new_execution_state_items]))
        new_releases_by_commit = [
            {
                version: max(
                    [
                        item
                        for item in new_execution_state_items
                        if item.version == version
                    ],
                    key=lambda x: x.commit.date_utc,
                )
            }
            for version in unique_new_releases
        ]
        # extract all commit values from new_releases_by_commit
        commits_pending = [
            list(item.values())[0].commit.sha for item in new_releases_by_commit
        ]

        # update execution state status as PENDING for commits that are the most recent release data,
        # and SKIPPED for commits that are not the most recent release data
        # ⚠️ see TODOs above regarding all commits from state being marked. This logic should process all commits 
        # that have null status, and mark them as PENDING or SKIPPED accordingly
        for item in new_execution_state_items:

            # mark item processing status
            if item.commit.sha in commits_pending:
                item.status = "PENDING"
            else:
                item.status = "SKIPPED"

            # set  execution_date_utc

        

    return


if __name__ == "__main__":
    lambda_handler({}, None)
