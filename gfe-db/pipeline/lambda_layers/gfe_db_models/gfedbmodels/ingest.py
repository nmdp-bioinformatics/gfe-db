import os
import logging
from typing import List, Dict, Union
import multiprocessing
from concurrent.futures import ThreadPoolExecutor, as_completed
from requests.exceptions import HTTPError
from time import sleep
from .types import (
    SourceConfig,
    RepositoryConfig,
    TargetMetadataConfig,
    Commit,
    ExecutionStateItem,
    ExecutionDetailsConfig,
)
from .utils import (
    cache_pickle,
    read_s3_json,
    sort_execution_state_items,
    filter_nulls,
    get_repo_asset,
    find_text,
)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


# TODO use method SourceConfig().from_s3(bucket, key) instead
def read_source_config(s3_client, bucket, key):
    data = read_s3_json(s3_client, bucket, key)
    return SourceConfig(**data)


def process_execution_state_item(
    timestamp: str,
    commit: Dict[str, str],
    repository_config: RepositoryConfig,
    target_metadata_config: TargetMetadataConfig,
    token: str = None,
    limit: int = None,
) -> Dict[str, str]:
    errors = 0
    sha = commit["sha"]

    for config in target_metadata_config.items:
        try:
            logger.info(
                f"Getting release version for sha {sha} from {config.asset_path}"
            )
            release_version = get_release_version_for_commit(
                commit=commit,
                owner=repository_config.owner,
                repo=repository_config.name,
                token=token,
                asset_path=config.asset_path,
                metadata_regex=config.metadata_regex,
            )
            logger.info(f"Found release version {release_version} ({sha})")

            result = (
                None,
                {
                    "created_utc": timestamp,
                    "repository": repository_config,
                    "commit": Commit(**commit),
                    "execution": ExecutionDetailsConfig(
                        version=release_version,
                        status="NOT_PROCESSED",
                        date_utc=None,
                        input_parameters=None,
                    ),
                },
            )

        except HTTPError as e:
            # This is because Allelelist.txt for certain commits doesn't contain the release version or name
            # Need to find another file that indicates the release version should be small
            errors += 1
            logger.info(f"Commit {sha} failed: {e} for {config.asset_path}")

            # Throw error if all possible asset paths have been tried
            if errors == len(target_metadata_config.items):
                logger.error(f"Max errors reached, this sha must be skipped: {sha}")

                # return the sha if all asset paths have been tried
                result = (sha, None)
                return result
            else:
                continue

        # return error count and increment outside this function

        # TODO deserialize to ExecutionStateItem, use as method
        return result


def parallel_process_execution_state_items(
    timestamp: str,
    commits: List[Dict[str, str]],
    repository_config: RepositoryConfig,
    target_metadata_config: TargetMetadataConfig,
    token: str = None,
    limit: int = None,
):
    execution_state_items = []
    num_cores = multiprocessing.cpu_count()
    num_threads = max(1, num_cores - 1)  # Reserve one core for other processes
    num_threads = min(6, num_cores)  # limit threads to avoid GitHub API rate limit

    # Create a ThreadPoolExecutor with the specified number of threads
    with ThreadPoolExecutor(max_workers=num_threads) as executor:

        # Submit the process_commit function for each commit to the executor
        futures = [
            executor.submit(
                process_execution_state_item,
                timestamp,
                commit,
                repository_config,
                target_metadata_config,
                token,
            )
            for commit in commits[:limit]
        ]

        # Collect the results as they complete
        results = [future.result() for future in as_completed(futures)]

        # Separate the results into execution state items and error shas
        execution_state_items = [
            result[1] for result in results if result[1] is not None
        ]
        error_shas = [result[0] for result in results if result[0] is not None]

    return (
        error_shas,
        [
            ExecutionStateItem(**item)
            for item in sort_execution_state_items(filter_nulls(execution_state_items))
        ],
    )


# limit is int or None
# @cache_pickle
def process_execution_state_items(
    timestamp: str,
    commits: List[Dict[str, str]],
    repository_config: RepositoryConfig,
    target_metadata_config: TargetMetadataConfig,
    token: str = None,
    limit: None = None,
    parallel: str = False,
) -> List[Dict[str, str]]:
    if parallel == True:
        if limit:
            logger.warning("'limit' will not work if parallel processing is enabled")
        return parallel_process_execution_state_items(
            timestamp=timestamp,
            commits=commits,
            repository_config=repository_config,
            target_metadata_config=target_metadata_config,
            token=token,
            limit=limit,
        )
    else:
        execution_state_items = []
        for commit in commits[:limit]:
            execution_state_items.append(
                process_execution_state_item(
                    timestamp=timestamp,
                    commit=commit,
                    repository_config=repository_config,
                    target_metadata_config=target_metadata_config,
                    token=token,
                    limit=limit,
                )
            )

        return [
            ExecutionStateItem(**item)
            for item in sort_execution_state_items(filter_nulls(execution_state_items))
        ]


def get_release_version_for_commit(
    commit: Union[Commit, dict],
    owner: str,
    repo: str,
    token: str,
    asset_path: str,
    metadata_regex: str,
) -> int:

    try:
        sha = commit["sha"]
    except:
        sha = commit.sha
    allele_list = get_repo_asset(
        owner=owner, repo=repo, token=token, path=asset_path, commit_sha=sha
    )

    release_version = find_text(metadata_regex, allele_list)

    if release_version is None:

        # TODO 2/12/24 save these shas (`fatal: reference is not a tree`) for debugging, need to get the file contents
        raise Exception(f"Release version not found for commit {sha}")

    return int(release_version.replace(".", "")[:4])


### debug ###
if __name__ == "__main__":
    import sys

    sys.path.append(os.environ["GFEDBMODELS_PATH"])
    from pathlib import Path
    import json
    import boto3
    from datetime import datetime
    from gfedbmodels.types import ExecutionStatus
    from gfedbmodels.utils import get_utc_now

    s3 = boto3.client("s3")

    utc_now = get_utc_now()

    GITHUB_REPOSITORY_OWNER = os.environ["GITHUB_REPOSITORY_OWNER"]
    GITHUB_REPOSITORY_NAME = os.environ["GITHUB_REPOSITORY_NAME"]
    GITHUB_PERSONAL_ACCESS_TOKEN = os.environ["GITHUB_PERSONAL_ACCESS_TOKEN"]

    test_path = json.loads(
        (
            Path(__file__).parent.parent.parent.parent
            / "functions"
            / "check_source_update"
            / "most-recent-commits.json"
        ).read_text()
    )
    commits_with_releases = []
    for commit in test_path:

        # Get data source configuration
        source_repo_config = read_source_config(
            s3_client=s3,
            bucket=os.environ["DATA_BUCKET_NAME"],
            key=os.environ["PIPELINE_SOURCE_CONFIG_S3_PATH"],
        ).repositories[f"{GITHUB_REPOSITORY_OWNER}/{GITHUB_REPOSITORY_NAME}"]

        # Loop through available file assets containing release version metadata
        for asset_config in source_repo_config.target_metadata_config.items:

            # Get the release version for the commit by examining file asset contents
            release_version = get_release_version_for_commit(
                commit=commit,
                owner=GITHUB_REPOSITORY_OWNER,
                repo=GITHUB_REPOSITORY_NAME,
                token=GITHUB_PERSONAL_ACCESS_TOKEN,
                asset_path=asset_config.asset_path,
                metadata_regex=asset_config.metadata_regex,
            )
            logger.info(
                f'Found release version {release_version} for commit {commit["sha"]}'
            )

            # Build the execution object to be stored in the state table (`execution__*` fields)
            execution_detail = ExecutionDetailsConfig(
                **{"version": release_version, "status": ExecutionStatus.PENDING}
            )

            # Build the repository object to be stored in the state table (`repository__*` fields)
            repository_config = RepositoryConfig(
                **{
                    "owner": GITHUB_REPOSITORY_OWNER,
                    "name": GITHUB_REPOSITORY_NAME,
                    "url": f"https://github.com/{GITHUB_REPOSITORY_OWNER}/{GITHUB_REPOSITORY_NAME}",
                }
            )

            # Assemble the execution state item for the new commit
            execution_state_item = ExecutionStateItem(
                created_utc=utc_now,
                execution=execution_detail,
                repository=repository_config,
                commit=Commit.from_response_json(commit),
            )
            commits_with_releases.append(execution_state_item)

            # break the loop if successful
            break
