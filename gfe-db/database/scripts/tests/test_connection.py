import os
import logging
import json
import boto3
from neo4j import GraphDatabase

logger = logging.getLogger()
logger.setLevel(logging.INFO)
STAGE = os.environ["STAGE"]
APP_NAME = os.environ["APP_NAME"]
AWS_REGION = os.environ["AWS_REGION"]

session = boto3.session.Session(region_name=AWS_REGION)
ssm = session.client("ssm")
secrets = session.client("secretsmanager")

# get neo4j uri and password from SSM Parameter Store
uri = ssm.get_parameter( Name=f"/{APP_NAME}/{STAGE}/{AWS_REGION}/Neo4jUri" )["Parameter"]["Value"]
# uri = ssm.get_parameter( Name=f"/{APP_NAME}/{STAGE}/{AWS_REGION}/Neo4jPrivateIp" )["Parameter"]["Value"]

# 'https://gfe-db.cloudftl.org:7473/browser/' => neo4j+s://gfe-db.cloudftl.org:7687 uri = "/".join(uri.replace("https://",
# "neo4j+s://").replace(":7473", ":7687").split("/")[:-2])
# uri = "/".join(uri.replace("https://", "neo4j://").split("/")) + ":7687"
# uri_scheme = "bolt" + "://"
# uri_port = ":" + "7687"
# uri = uri_scheme + uri + uri_port
logger.info(f"Neo4j URI: {uri}")

# /gfe-db/dev/us-east-1/Neo4jCredentialsSecretArn
auth_arn = ssm.get_parameter( Name=f"/{APP_NAME}/{STAGE}/{AWS_REGION}/Neo4jCredentialsSecretArn" )["Parameter"]["Value"]

# get secret from arn
auth = json.loads(secrets.get_secret_value(SecretId=auth_arn)["SecretString"])

if __name__ == "__main__":
        graphdb = GraphDatabase.driver(uri, auth=(auth["NEO4J_USERNAME"], auth["NEO4J_PASSWORD"]))
        print(graphdb)

        with graphdb as driver:
            records, _, _ = driver.execute_query(f'SHOW CONSTRAINTS;', database_="neo4j")

        print(records)