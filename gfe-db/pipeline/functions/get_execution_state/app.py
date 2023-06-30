"""
Checks a GitHub repository for new commits and triggers data ingestion. This function processes
only the releases that it finds. To process specific releases, use a different method.

Note: this function is only responsible for checking and processing the most recent commits. It is not responsible for 
syncing state. If old items are deleted on the Execution state table while the most recent commits remain, 
this function will not reprocess the deleted items.
"""
import os

if __name__ != "app":
    import sys

    # for dev, local path to gfe-db modules
    # ./gfe-db/pipeline/lambda_layers/gfe_db_models (use absolute path)
    sys.path.append(os.environ["GFEDBMODELS_PATH"])
import logging
import json

logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Environment
APP_NAME = os.environ["APP_NAME"]
STAGE = os.environ["STAGE"]


def lambda_handler(event, context):
    logger.info(json.dumps(event))

    # TODO validate input
    # TODO Get state for commit in input
    # TODO Return state
    # TODO return the SQS message receipt in case it needs to be returned to the queue if the state machine fails

    return


if __name__ == "__main__":
    from pathlib import Path

    event = json.loads((Path(__file__).parent / "event.json").read_text())

    lambda_handler(event, None)
