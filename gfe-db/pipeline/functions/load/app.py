
import logging
import json
import boto3

logger = logging.getLogger()
logger.setLevel(logging.INFO)

def lambda_handler(event, context):
    """Invoke SSM Run Command for server side loading on Neo4j"""

    logger.info(json.dumps(event))

    return


if __name__ == "__main__":

    path = '/Users/ammon/Documents/00-Projects/nmdp-bioinformatics/02-Repositories/gfe-db/gfe-db/pipeline/functions/trigger/src/event.json'
    # path = '/Users/ammon/Documents/00-Projects/nmdp-bioinformatics/02-Repositories/gfe-db/gfe-db/pipeline/functions/trigger/src/schedule-event.json'

    with open(path, "r") as file:
        event = json.load(file)

    lambda_handler(event,"")
