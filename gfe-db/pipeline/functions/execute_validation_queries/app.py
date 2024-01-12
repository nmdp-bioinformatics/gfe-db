"""
This function executes validation queries against the Neo4j database and returns the results.
If USE_PRIVATE_SUBNET is true, this function will run inside a VPC and private subnet. 
If USE_PRIVATE_SUBNET is false, this function will run outside a VPC and in a public subnet.
"""
import os
import logging
import json
import boto3
from neo4j import GraphDatabase

logger = logging.getLogger()
logger.setLevel(logging.INFO)

AWS_REGION = os.environ["AWS_REGION"]
STAGE = os.environ["STAGE"]
APP_NAME = os.environ["APP_NAME"]

# get neo4j uri and password from SSM Parameter Store
session = boto3.session.Session(region_name=AWS_REGION)
ssm = session.client("ssm")
secrets = session.client("secretsmanager")

uri = ssm.get_parameter(
    Name=f"/{APP_NAME}/{STAGE}/{AWS_REGION}/Neo4jUri"
)["Parameter"]["Value"]
logger.info(f"Found uri: {uri}")

auth = json.loads(secrets.get_secret_value(SecretId=f'/{APP_NAME}/{STAGE}/{AWS_REGION}/Neo4jCredentials')["SecretString"])
logger.info(f'Found auth: {auth["NEO4J_USERNAME"]}')

logger.info(f"Connecting to Neo4j at {uri}")
graphdb = GraphDatabase.driver(uri, auth=(auth["NEO4J_USERNAME"], auth["NEO4J_PASSWORD"]))

# TODO get database name from SSM Parameter Store
def lambda_handler(event, context):

    logger.info(json.dumps(event))

    with graphdb as driver:

        # node counts
        node_counts = []
        for node in nodes:
            records, _, _ = driver.execute_query(f'MATCH (n:{node}) RETURN count(n) as count;', database_="gfedb")
            node_counts.append({
                "node": node,
                "count": records[0].data()['count']
            })

        # HAS_IPD_ALLELE relationship releases property release counts
        has_ipd_allele_release_counts = execute_query(driver, has_ipd_allele_release_counts_cql)

        # IPD_Accession node release counts
        ipd_accession_release_counts = execute_query(driver, ipd_accession_release_counts_cql)

    payload = {
        "node_counts": node_counts,
        "has_ipd_allele_release_counts": has_ipd_allele_release_counts,
        "ipd_accession_release_counts": ipd_accession_release_counts
    }
    return payload

nodes = [
    "GFE",
    "IPD_Accession",
    "IPD_Allele",
    "Sequence",
    "Feature",
    "Submitter",
]

has_ipd_allele_release_counts_cql = """MATCH (:GFE)-[r:HAS_IPD_ALLELE]->(:IPD_Allele)
WITH r, apoc.coll.toSet(r.releases) as releases
UNWIND releases as release_version
RETURN DISTINCT release_version, count(release_version) as count
ORDER BY release_version;"""

ipd_accession_release_counts_cql = """MATCH ()-[r:HAS_IPD_ACCESSION]->() RETURN DISTINCT r.release as release_version, count(r.release) as count;"""

# # too slow
# node_counts_cql = """MATCH (gfe:GFE)
# MATCH (ipd:IPD_Accession)
# MATCH (ipda:IPD_Allele)
# MATCH (seq:Sequence)
# MATCH (f:Feature)
# MATCH (sub:Submitter)
# RETURN 
#     count(gfe),
#     count(ipd),
#     count(ipda),
#     count(seq),
#     count(f),
#     count(sub);"""

def execute_query(driver, query, database="gfedb"):
    records, _, _ = driver.execute_query(query, database_=database)
    return [record.data() for record in records]

if __name__ == "__main__":
    from pathlib import Path

    event_path = Path(__file__).parent / "post-execution-event.json"

    with open(event_path, "r") as file:
        event = json.load(file)

    lambda_handler(event,"")
