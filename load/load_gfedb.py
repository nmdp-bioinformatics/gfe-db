import os
import logging
import time
import base64
import json
import requests
import boto3

logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Environment variables
root = os.environ["ROOT"]
neo4j_dir = os.environ["NEO4J_DIR"]
cypher_dir = "/".join([neo4j_dir, "cypher"])
load_script = os.environ["SCRIPT"]
s3_bucket = os.environ["GFE_BUCKET"]
release = os.environ["RELEASES"]

host = os.environ["NEO4J_HOST"]
username = os.environ["NEO4J_USERNAME"]
password = os.environ["NEO4J_PASSWORD"]
protocol = 'http'
port = "7474"
endpoint = "db/neo4j/tx/commit"
url = f'{protocol}://{host}:{port}/{endpoint}'

s3_urls = [
    's3://gfe-db-4498/data/RELEASE/csv/all_groups.RELEASE.csv',
    's3://gfe-db-4498/data/RELEASE/csv/all_cds.RELEASE.csv',
    's3://gfe-db-4498/data/RELEASE/csv/all_features.RELEASE.csv',
    's3://gfe-db-4498/data/RELEASE/csv/gfe_sequences.RELEASE.csv',
    's3://gfe-db-4498/data/RELEASE/csv/all_alignments.RELEASE.csv'
]


def generate_presigned_urls(s3_urls, expire=3600):
    """Accepts a list of S3 URLs or paths and returns
    a dictionary of pre-signed URLs for each"""
    
    logger.info(f"Generating pre-signed URLs...")
    
    s3_urls = [s3_urls] if not isinstance(s3_urls, list) else s3_urls
    
    presigned_urls = {}
    
    for s3_url in s3_urls:
        
        i = 2 if "s3://" in s3_url else 0
        
        bucket = s3_url.split("/")[i]
        key = "/".join(s3_url.split("/")[i + 1:])
        
        # Generate the URL to get 'key-name' from 'bucket-name'
        url = s3.generate_presigned_url(
            ClientMethod='get_object',
            Params={
                'Bucket': bucket,
                'Key': key
            },
            ExpiresIn=expire
        )
        
        presigned_urls[s3_url] = url
        
    return presigned_urls


def update_cypher(cypher_path):
    """Replaces instances of "file:///{csv_prefix}.RELEASE.csv" with
    an S3 pre-sign URL"""

    with open(cypher_path, "r") as file:
        cypher_script = file.read()

    for s3_url in s3_urls:

        csv_prefix = s3_url.split("/")[-1].split(".")[0]
        cypher_script = cypher_script.replace(f'file:///{csv_prefix}.RELEASE.csv', presigned_urls[s3_url])
        
    return cypher_script


def run_cypher(cypher):
    
    payload = {
        "statements": [
            {
                "statement": cypher
            }
        ]
    }
    
    logger.debug(payload)
    
    # Headers
    headers = { 
        "Accept": "application/json;charset=UTF-8",
        "Content-Type": "application/json",
        "Authorization": f"Basic {base64.b64encode(':'.join([username, password]).encode()).decode()}"
    }

    # Send requests
    response = requests.post(
        url, 
        data=json.dumps(payload), 
        headers=headers)
    
    try:
        response_dict = json.loads(response.content)
      
        if len(response_dict['errors']) > 0:
            logger.error(response_dict)
        else:
            logger.info(response_dict)
        
        return response_dict

    except Exception as err:
        logger.error(f"Failed to load response from Neo4j server: {response.status_code}")
        #raise err
        return


if __name__ == "__main__":

    s3_urls = list(map(lambda x: x.replace("RELEASE", release), s3_urls))

    # Get the service client.
    s3 = boto3.client('s3')
    presigned_urls = generate_presigned_urls(s3_urls)

    cypher_path = "/".join([f'{cypher_dir}/{load_script}'])
    cypher_script = update_cypher(cypher_path)
    cypher = list(filter(lambda x: x != "\n", cypher_script.split(";")))
    cypher = list(map(lambda x: "".join([x, ";"]), cypher))

    #limit = 1
    start = time.time()
    for idx, statement in enumerate(cypher):
        print(f'Executing statement: {statement}')
        statement_start = time.time()
        response = run_cypher(statement)
        statement_end = time.time()
        statement_elapsed_time = statement_end - statement_start
        logger.info(f'Statement: {statement}\nTime elapsed: {statement_elapsed_time}\nResponse: {response}')
        print(f'Time elapsed: {statement_elapsed_time}\nResponse: {response}\n\n')
        #if limit and idx + 1 == limit:
        #    break
            
    end = time.time()
    time_elapsed = end - start
    print(f"Time elapsed: {time_elapsed} seconds")