import os
import logging
import json
import boto3

logger = logging.getLogger()
logger.setLevel(logging.INFO)

def lambda_handler(event, context):
    logger.info(json.dumps(event, indent=4))
    
    return


# Needed to serialize datetime objects in JSON responses
class DatetimeEncoder(json.JSONEncoder):
    """
    Helps convert datetime objects to pure strings in AWS service API responses. Does not
    convert timezone information.

    Extend `json.JSONEncoder`. 
    """

    def default(self, obj):
        try:
            return super().default(obj)
        except TypeError:
            return str(obj)



if __name__ == "__main__":

    path = '/Users/ammon/Documents/00-Projects/nmdp-bioinformatics/02-Repositories/gfe-db/gfe-db/pipeline/functions/invoke_load_script/event.json'

    with open(path, "r") as file:
        event = json.load(file)

    lambda_handler(event,"")
