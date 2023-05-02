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
    InputParameters,
    ExecutionStateItem,
    RepositoryConfig,
    SourceConfig,
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
GITHUB_REPOSITORY_OWNER = "ANHIG"  # os.environ["GITHUB_REPOSITORY_OWNER"]
GITHUB_REPOSITORY_NAME = "IMGTHLA"  # os.environ["GITHUB_REPOSITORY_NAME"]
DATA_BUCKET_NAME = os.environ["DATA_BUCKET_NAME"]
PIPELINE_CONFIG_S3_PATH = os.environ["PIPELINE_CONFIG_S3_PATH"]


if __name__ == "__main__":
    # Paths
    # TODO arg
    output_dir = Path(f"{APP_NAME}/pipeline/config")

    # Fetch all commits from repo using GitHub API, will be cached
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
    max_errors = 5
    limit = None
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

    # Excluding old commits with no useful information
    excluded_commit_shas = [
        "08e0ef9f5c6aade40df681821a0b9caef439fe3a",
        "6ad21b61dee3689c5ae68370d635c5ede483c851",
        "79d13ceb388eb9dacc9e166be18cce9373f7fd1d",
        "9f35f8fe8a2e25bb076e588e65389cac16a8ed2f",
        "785c913f2d42abd68bcdf630ce2f58ee9b9c2579",
        "efc06e88b56d1e6e44661ec45f192dc1186a30ad",
    ]
    commits = [commit for commit in commits if commit["sha"] not in excluded_commit_shas]

    # Build ExecutionStateItem list using thread pooling    
    execution_state_items = process_execution_state_items(
        commits=commits,
        asset_configs=asset_configs,
        limit=limit,
        parallel=True,
    )

    # Build RepositoryConfig
    repository_path = f"{GITHUB_REPOSITORY_OWNER}/{GITHUB_REPOSITORY_NAME}"
    tracked_assets = ["hla.dat", "msf/"]

    source_config = SourceConfig(
        **{
            "created_at_utc": utc_now,
            "updated_at_utc": utc_now,
            "repositories": {
                repository_path: RepositoryConfig(
                    **{
                        "owner": GITHUB_REPOSITORY_OWNER,
                        "name": GITHUB_REPOSITORY_NAME,
                        "url": f"https://github.com/{repository_path}",
                        "tracked_assets": tracked_assets,
                        # TODO fetch default input parameters from S3 so they are easy to update
                        "default_input_parameters": InputParameters(
                            **{
                                "align": os.environ.get("ALIGN", False),
                                "kir": os.environ.get("KIR", False),
                                "mem_profile": os.environ.get("MEM_PROFILE", False),
                                "limit": os.environ.get("LIMIT", None),
                            }
                        ),
                        "execution_state": [
                            ExecutionStateItem(**item)
                            for item in execution_state_items
                        ],
                    }
                )
            },
        }
    )

    # # convert utc_now string to datetime object using YYYYMMDDHHMM format
    # utc_now_version = datetime.strptime(utc_now, "%Y-%m-%dT%H:%M:%SZ").strftime(
    #     "%Y%m%d-%H%M"
    # )

    # TODO get from argparse
    # write SourceConfig locally
    with open(output_dir / f"source-config.json", "w") as f:
        json.dump(source_config.dict(), f, indent=4)

exit(0)
