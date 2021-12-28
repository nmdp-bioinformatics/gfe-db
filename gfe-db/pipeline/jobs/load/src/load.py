import os
import logging
import time
import base64
import ast
import json
import requests
import boto3

# TODO: update headers to include "X-Stream": "true"

logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Environment variables
root = os.environ["ROOT"]
cypher_dir = os.environ["CYPHER_DIR"]
load_script = os.environ["LOAD_SCRIPT"]
s3_bucket = os.environ["GFE_BUCKET"]
release = os.environ["RELEASES"]
alignments = ast.literal_eval(os.environ["ALIGN"])

# TODO: add KIR
kir = ast.literal_eval(os.environ["KIR"])

host = os.environ["NEO4J_HOST"]
username = os.environ["NEO4J_USERNAME"]
password = os.environ["NEO4J_PASSWORD"]
protocol = 'http'
port = "7474"
endpoint = "db/neo4j/tx/commit"
url = f'{protocol}://{host}:{port}/{endpoint}'

s3_urls = [
    f's3://{s3_bucket}/data/{release}/csv/all_groups.{release}.csv',
    f's3://{s3_bucket}/data/{release}/csv/all_cds.{release}.csv',
    f's3://{s3_bucket}/data/{release}/csv/all_features.{release}.csv',
    f's3://{s3_bucket}/data/{release}/csv/gfe_sequences.{release}.csv'
]

if alignments:
    s3_urls.append(f's3://{s3_bucket}/data/{release}/csv/all_alignments.{release}.csv')

if kir: 
    s3_urls.append(f's3://{s3_bucket}/data/{release}/csv/all_kir.{release}.csv')

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


def clean_cypher(cypher_script):
    """Cleans Cypher script to remove comments, new lines, add semi-colons etc."""

    cypher = list(filter(lambda x: x != "\n", cypher_script.split(";")))
    cypher = list(map(lambda x: "".join([x, ";"]), cypher))
    cypher = list(map(lambda x: x.replace("\n", " ").strip(), cypher))
    cypher = list(filter(lambda x: "// RETURN" not in x, cypher))

    if not alignments:
        cypher = list(filter(lambda x: "all_alignments" not in x, cypher))

    return cypher

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
    logger.info(f'Running query:\n{cypher}')
    response = requests.post(
        url, 
        data=json.dumps(payload), 
        headers=headers)
    
    try:
        response_dict = json.loads(response.content)
      
        if len(response_dict['errors']) > 0:
            for error in response_dict['errors']:
                logger.error(f'Error loading Neo4j: {error}')
        else:
            logger.info(response_dict)
        
        return response_dict

    except Exception as err:
        logger.error(f"Failed to load response from Neo4j server: {response.status_code}")
        raise err


if __name__ == "__main__":

    # s3_urls = list(map(lambda x: x.replace("RELEASE", release), s3_urls))

    s3 = boto3.client('s3')

    # Not working in us-east-2, possibly because of endpoint_url formatting
    # See: https://github.com/boto/boto3/issues/1982
    presigned_urls = generate_presigned_urls(s3_urls)
    time.sleep(10)

    cypher_path = "/".join([f'{cypher_dir}/{load_script}'])
    cypher_script = update_cypher(cypher_path)
    cypher = clean_cypher(cypher_script)
    logger.info(f'Cypher script: {cypher}')

    start = time.time()
    for idx, statement in enumerate(cypher):
        logger.info(f'Executing statement: {statement}')
        statement_start = time.time()
        response = run_cypher(statement)
        statement_end = time.time()
        statement_elapsed_time = round(statement_end - statement_start, 2)
        logger.info(f'Loaded in {statement_elapsed_time} s\nResponse: {response}\n')
            
    end = time.time()
    time_elapsed = round(end - start, 2)
    logger.info(f"Total time elapsed: {time_elapsed} seconds")