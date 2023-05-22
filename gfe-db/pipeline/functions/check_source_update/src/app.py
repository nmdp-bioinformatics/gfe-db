"""
Checks a GitHub repository for new commits and triggers data ingestion. This function processes
only the releases that it finds. To process specific releases, use a different method.
"""
import os
import logging
from decimal import Decimal
from datetime import datetime, timedelta
utc_now = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
import json
import boto3
from utils.constants import (
    GITHUB_REPOSITORY_OWNER,
    GITHUB_REPOSITORY_NAME,
    execution_state_table_fields
)
from utils.types import (
    str_to_datetime,
    str_from_datetime,
    ExecutionStateItem, 
    RepositoryConfig,
    Commit,
    ExecutionDetailsConfig,
    ExecutionPayloadItem
)
from utils.utils import (
    read_source_config,
    restore_nested_json, 
    list_commits, 
    get_release_version_for_commit,
    flatten_json_records,
    flatten_json,

    cache_json,
    cache_pickle
)

# set up logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Environment
AWS_REGION = os.environ["AWS_REGION"]
APP_NAME = os.environ["APP_NAME"]
STAGE = os.environ["STAGE"]
PIPELINE_SOURCE_CONFIG_S3_PATH = os.environ["PIPELINE_SOURCE_CONFIG_S3_PATH"]
GITHUB_REPOSITORY_OWNER = os.environ["GITHUB_REPOSITORY_OWNER"]
GITHUB_REPOSITORY_NAME = os.environ["GITHUB_REPOSITORY_NAME"]

# AWS clients
session = boto3.Session(region_name=AWS_REGION)
ssm = session.client("ssm")
dynamodb = session.resource("dynamodb", region_name=AWS_REGION)

# Get SSM Parameters
table_name = ssm.get_parameter(
    Name=f'/{os.environ["APP_NAME"]}/{os.environ["STAGE"]}/{os.environ["AWS_REGION"]}/ExecutionStateTableName'
)["Parameter"]["Value"]

state_machine_arn = ssm.get_parameter(
    Name=f'/{os.environ["APP_NAME"]}/{os.environ["STAGE"]}/{os.environ["AWS_REGION"]}/UpdatePipelineArn'
)["Parameter"]["Value"]

data_bucket_name = ssm.get_parameter(
    Name=f'/{os.environ["APP_NAME"]}/{os.environ["STAGE"]}/{os.environ["AWS_REGION"]}/DataBucketName'
)["Parameter"]["Value"]

# Get data source configuration
source_repo_config = (
    read_source_config(data_bucket_name, PIPELINE_SOURCE_CONFIG_S3_PATH)
    .repositories[f"{GITHUB_REPOSITORY_OWNER}/{GITHUB_REPOSITORY_NAME}"]
)


def lambda_handler(event, context):
    logger.info(json.dumps(event))

    # Get items from state table
    table = dynamodb.Table(table_name)

    execution_state = get_execution_state(table)

    # # ⬇ TESTING uncomment before deploying ⬇
    # test_items_to_delete = execution_state[:5]

    # # # delete from table using commit.sha
    # # # ValidationError
    # # for item in test_items_to_delete:
    # #     # table.delete_item(Key={"commit.sha": item.commit.sha})

    # del execution_state[:5]
    # # ⬆ TESTING uncomment before deploying ⬆

    # Get the most recent commits from github since the most recent commit date retrieved from DynamoDB
    commits = get_most_recent_commits(execution_state)

    if not commits:
        logger.info("No new commits found")
        return

    # TODO Build a list of ExecutionStateItems using defaults for execution, repository and adding the commit BOOKMARK (1)
    commits_with_releases = []
    for commit in commits:
        sha = commit["sha"]
        for asset_config in source_repo_config.target_metadata_config.items:
            try:
                release_version = get_release_version_for_commit(commit, **asset_config.dict())
                execution_detail = ExecutionDetailsConfig(**{"version": release_version, "status": "NOT_PROCESSED"})
                repository_config = RepositoryConfig(**{
                    "owner": GITHUB_REPOSITORY_OWNER, 
                    "name": GITHUB_REPOSITORY_NAME,
                    "url": f"https://github.com/{GITHUB_REPOSITORY_OWNER}/{GITHUB_REPOSITORY_NAME}",
                    "default_input_parameters": source_repo_config.default_input_parameters,

                })
                commit = Commit.from_response_json(commit)
                execution_state_item = ExecutionStateItem(
                    updated_utc=utc_now,
                    execution=execution_detail,
                    repository=repository_config,
                    commit=commit
                )
                commits_with_releases.append(execution_state_item)
                # break the loop if successful
                break
            except Exception as e:
                logger.info(f"Error getting release version for commit {sha}: {e}")

    ### Trigger the build process for each release with the most recent commit for that version ###
    # 1) Mark the most recent commit for each release as PENDING
    # TODO add input parameters, take into account that they will come from either the default or from user input through event
    input_parameters = source_repo_config.default_input_parameters
    pending_commits = [ 
        update_execution_state_item(commit, status="PENDING", input_parameters=input_parameters) \
            for item in select_most_recent_commit_for_release(commits_with_releases) \
                for commit in item.values()
    ]


    # 2) Mark the older commits for each release as SKIPPED
    skipped_commits = [ update_execution_state_item(commit, "SKIPPED") for commit in commits_with_releases if commit not in pending_commits ]

    # 3) Update the state table with the new commits, order by commit.date_utc descending
    new_execution_state = pending_commits + skipped_commits
    new_execution_state = sorted(new_execution_state, key=lambda x: x.commit.date_utc, reverse=False)

    # 4) Flatten and load the processed records to the state table
    items = [ 
        flatten_json(
            data=item.dict(),
            select_fields=execution_state_table_fields) \
                for item in new_execution_state
        ]
    # for item in new_execution_state:

    #     # flatten the item and select only the fields in the table
    #     items.append(flatten_json(
    #         data=item.dict(),
    #         select_fields=execution_state_table_fields))

    with table.batch_writer() as batch:
        logger.info(f"Loading {len(items)} items to {table_name}")
        for item in items:
            batch.put_item(Item=item)

    # 5) Return pending commits to the state machine for further processing
    # TODO add input parameters
    execution_payload = [ ExecutionPayloadItem.from_execution_state_item(item).dict() for item in pending_commits ] 
    return execution_payload

@cache_pickle
def get_execution_state(table):
    # Retrieve execution state from table
    items = table.scan()["Items"]
    # convert Decimal types to int
    from decimal import Decimal
    items = [{k: int(v) if isinstance(v, Decimal) else v for k, v in item.items()} for item in items]
    items = sorted(items, key=lambda x: x["commit.date_utc"], reverse=True)

    # TODO Deserialize and repack the items
    return [ ExecutionStateItem(**restore_nested_json(item)) for item in items ]

@cache_json
def get_most_recent_commits(execution_state):
    # 1) Get the most recent commit date from DynamoDB using max(), add one second to it so the same commit is not returned
    last_commit_date = max([ str_to_datetime(item.commit.date_utc) for item in execution_state ])
    
    # add minor offset to avoid duplicate commits
    since = str_from_datetime(last_commit_date + timedelta(seconds=1))

    # 2) Get the most recent commits from GitHub using since=<date> parameter
    return list_commits(GITHUB_REPOSITORY_OWNER, GITHUB_REPOSITORY_NAME, since=since)

def select_most_recent_commit_for_release(commits: list[ExecutionStateItem]):
    # group by release version and get most recent by commit date (max date_utc)
    unique_new_releases = list(set([item.execution.version for item in commits]))
    return [
        {
            version: max(
                [
                    item
                    for item in commits
                    if item.execution.version == version
                ],
                key=lambda x: x.commit.date_utc,
            )
        }
        for version in unique_new_releases
    ]

def update_execution_state_item(execution_state_item: ExecutionStateItem, status: str, input_parameters: dict = None):
    execution_state_item.execution.status = status

    if input_parameters is not None and status == "PENDING":
        execution_state_item.execution.input_parameters = input_parameters

    execution_state_item.execution.date_utc = str_from_datetime(datetime.utcnow())
    return execution_state_item



if __name__ == "__main__":
    lambda_handler({}, None)
