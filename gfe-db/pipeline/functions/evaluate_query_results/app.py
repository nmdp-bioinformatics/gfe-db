"""
Query result evaluation:
For post-load invocations the script appends additional metadata to the payload to indicate whether the load was successful.
The success condition is defined by the following criteria:
- The number of nodes in the database is greater than the number of nodes in the pre-load invocation.
- The post-load release version for the State Machine execution matches the release version in the database.
- The post-load number of unique release versions is greater than the number of unique release versions in the pre-load invocation.
- If a limit is specified, the number of GFE nodes for the specific release matches the limit
"""
import os
import logging
import json

logger = logging.getLogger()
logger.setLevel(logging.INFO)

AWS_REGION = os.environ["AWS_REGION"]
STAGE = os.environ["STAGE"]
APP_NAME = os.environ["APP_NAME"]

def lambda_handler(event, context):

    release_version = event['input']['version']
    query_results = event['validations']['queries']

    # Release has been added to the database
    unique_releases_in_db_pre_load = sorted([ int(item['release_version']) for item in query_results['pre']['has_ipd_allele_release_counts'] ])
    unique_releases_in_db_post_load = sorted([ int(item['release_version']) for item in query_results['post']['has_ipd_allele_release_counts'] ])
    is_release_version_loaded = set(unique_releases_in_db_post_load) - set(unique_releases_in_db_pre_load) == set([int(release_version)])

    # Number of nodes in the database has increased
    node_counts_pre_load = sum(sorted([ item['count'] for item in query_results['pre']['node_counts'] ]))
    node_counts_post_load = sum(sorted([ item['count'] for item in query_results['post']['node_counts'] ]))
    have_node_counts_increased = node_counts_post_load > node_counts_pre_load

    # Number of unique release versions in the database has increased by one
    num_unique_releases_in_db_post_load = len(unique_releases_in_db_post_load)
    num_unique_releases_in_db_pre_load = len(unique_releases_in_db_pre_load)
    has_unique_release_count_increased_by_1 = num_unique_releases_in_db_post_load == num_unique_releases_in_db_pre_load + 1

    is_load_successful = (
        is_release_version_loaded
        and have_node_counts_increased
        and has_unique_release_count_increased_by_1
    )

    payload = {
        "is_load_successful": {
            "value": is_load_successful,
            "details": {
                "is_release_version_loaded": {
                    "value": is_release_version_loaded,
                    "details": {
                        "unique_releases_in_db_pre_load": unique_releases_in_db_pre_load,
                        "unique_releases_in_db_post_load": unique_releases_in_db_post_load
                    }
                },
                "have_node_counts_increased": {
                    "value": have_node_counts_increased,
                    "details": {
                        "node_counts_pre_load": node_counts_pre_load,
                        "node_counts_post_load": node_counts_post_load
                    },
                },
                "has_unique_release_count_increased_by_1": {
                    "value": has_unique_release_count_increased_by_1,
                    "details": {
                        "num_unique_releases_in_db_pre_load": num_unique_releases_in_db_pre_load,
                        "num_unique_releases_in_db_post_load": num_unique_releases_in_db_post_load
                    }
                }
            }
        }
    }

    logger.info(json.dumps(payload))

    return payload

if __name__ == "__main__":
    from pathlib import Path

    event_path = Path(__file__).parent / "false-event.json"

    with open(event_path, "r") as file:
        event = json.load(file)

    lambda_handler(event,"")
