import logging
from typing import List, Dict, Union
import multiprocessing
from concurrent.futures import ThreadPoolExecutor, as_completed
from .constants import (
    pipeline
)
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
    find_text
)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def read_source_config(s3_client, bucket, key):
    data = read_s3_json(s3_client, bucket, key)
    return SourceConfig(**data)

# def write_source_config(bucket, key, source_config: SourceConfig):
#     write_s3_json(bucket, key, source_config.dict())

def process_execution_state_item(
    timestamp: str,
    commit: Dict[str, str],
    repository_config: RepositoryConfig,
    target_metadata_config: TargetMetadataConfig,
    limit: int = None,
) -> Dict[str, str]:
    errors = 0
    sha = commit["sha"]

    for config in target_metadata_config.items:
        try:
            logger.info(
                f"Getting release version for sha {sha} from {config.asset_path}"
            )
            release_version = get_release_version_for_commit(commit, **config.dict())
            logger.info(f"Found release version {release_version} ({sha})")

            result = {
                "created_utc": timestamp,
                "repository": repository_config,
                "commit": Commit(**commit),
                "execution": ExecutionDetailsConfig(
                    version=release_version,
                    status="NOT_PROCESSED",
                    date_utc=None,
                    input_parameters=None,
                ),
            }

        except Exception as e:
            # This is because Allelelist.txt for certain commits doesn't contain the release version or name
            # Need to find another file that indicates the release version should be small
            errors += 1
            logger.info(f"Commit {sha} failed: {e}")

            # Throw error if all possible asset paths have been tried
            if errors == len(target_metadata_config.items):
                # logger.error(f"Max errors reached. Exiting loop.")
                raise e
            else:
                continue

        # return error count and increment outside this function

        # TODO deserialize to ExecutionStateItem
        return result
    
def parallel_process_execution_state_items(
    timestamp: str,
    commits: List[Dict[str, str]],
    repository_config: RepositoryConfig,
    target_metadata_config: TargetMetadataConfig,
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
            )
            for commit in commits[:limit]
        ]

        # Collect the results as they complete
        execution_state_items = [future.result() for future in as_completed(futures)]

    return [
        ExecutionStateItem(**item)
        for item in sort_execution_state_items(filter_nulls(execution_state_items))
    ]


# limit is int or None
@cache_pickle
def process_execution_state_items(
    timestamp: str,
    commits: List[Dict[str, str]],
    repository_config: RepositoryConfig,
    target_metadata_config: TargetMetadataConfig,
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
            limit=limit,
        )
    else:
        execution_state_items = []
        for commit in commits[:limit]:
            execution_state_items.append(
                process_execution_state_item(
                    commit=commit,
                    repository_config=repository_config,
                    target_metadata_config=target_metadata_config,
                    limit=limit,
                )
            )

        return [
            ExecutionStateItem(**item)
            for item in sort_execution_state_items(filter_nulls(execution_state_items))
        ]
    

def get_release_version_for_commit(commit: Union[Commit, dict], **kwargs) -> int:
    try:
        sha = commit["sha"]
    except:
        sha = commit.sha
    asset_path = kwargs["asset_path"]
    release_version_regex = kwargs["metadata_regex"]
    allele_list = get_repo_asset(
        pipeline.params.GitHubSourceRepository["owner"], pipeline.params.GitHubSourceRepository["name"], asset_path, sha
    )
    release_version = find_text(release_version_regex, allele_list)
    if release_version is None:
        raise Exception(f"Release version not found for commit {sha}")
    return int(release_version.replace(".", "")[:4])