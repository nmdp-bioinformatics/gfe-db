import os
import logging
import json
import boto3
from neo4j import GraphDatabase
from utils import format_uri

logger = logging.getLogger()
logger.setLevel(logging.INFO)

STAGE = os.environ["STAGE"]
APP_NAME = os.environ["APP_NAME"]
AWS_REGION = os.environ["AWS_REGION"]

# get neo4j uri and password from SSM Parameter Store
session = boto3.session.Session(region_name=AWS_REGION)
ssm = session.client("ssm")
secrets = session.client("secretsmanager")

# /APP_NAME/STAGE/AWS_REGION/Neo4jDatabaseEndpoint
uri = ssm.get_parameter(
    Name=f"/{APP_NAME}/{STAGE}/{AWS_REGION}/Neo4jUri"
)["Parameter"]["Value"]
# 'https://gfe-db.cloudftl.org:7473/browser/' => neo4j+s://gfe-db.cloudftl.org:7687
uri = format_uri(uri)

# /gfe-db/dev/us-east-1/Neo4jCredentialsSecretArn
auth_arn = ssm.get_parameter(
    Name=f"/{APP_NAME}/{STAGE}/{AWS_REGION}/Neo4jCredentialsSecretArn"
)["Parameter"]["Value"]

# get secret from arn
auth = json.loads(secrets.get_secret_value(SecretId=auth_arn)["SecretString"])

graphdb = GraphDatabase.driver(uri, auth=(auth["NEO4J_USERNAME"], auth["NEO4J_PASSWORD"]))

def lambda_handler(event, context):

    logger.info(json.dumps(event))

    # TODO write validation query
    # # Group and count nodes and edges for each release
    cql = "MATCH (n) RETURN count(n);"

    with graphdb as driver:
        records, _, _ = driver.execute_query(cql, database_="neo4j")

    return


if __name__ == "__main__":
    from pathlib import Path

    event_path = Path(__file__).parent / "event.json"

    with open(event_path, "r") as file:
        event = json.load(file)

    lambda_handler(event,"")
