#!/bin/bash

REGION=$(curl --silent http://169.254.169.254/latest/dynamic/instance-identity/document | jq -r .region)

# Set absolute paths
NEO4J_CYPHER_PATH_S3=config/neo4j/cypher
NEO4J_CYPHER_PATH=/var/lib/neo4j/cypher
NEO4J_IMPORT_PATH=/var/lib/neo4j/import

# Check for release argument
RELEASE=$1

if [[ -z $RELEASE ]]; then
    echo "Release version not found"
    exit 1
else
    echo "Starting load process for $RELEASE"
fi

# Get Neo4j Credentials
NEO4J_CREDENTIALS=$(aws secretsmanager get-secret-value \
    --region $REGION \
    --secret-id gfe-db-dev-Neo4jCredentials | jq -r '.SecretString')
NEO4J_USERNAME=$(echo $NEO4J_CREDENTIALS | jq -r '.NEO4J_USERNAME')
NEO4J_PASSWORD=$(echo $NEO4J_CREDENTIALS | jq -r '.NEO4J_PASSWORD')

# Get data bucket name
DATA_BUCKET_NAME=$(aws ssm get-parameters \
    --region $(curl -s http://169.254.169.254/latest/meta-data/placement/availability-zone | sed 's/\(.*\)[a-z]/\1/') \
    --names "/gfe-db/dev/us-east-1/DataBucketName" \
    | jq -r '.Parameters | map(select(.Version == 1))[0].Value')

if [[ -z $DATA_BUCKET_NAME ]]; then
    echo "S3 bucket not found."
    exit 1
else
    echo "Found S3 bucket: $DATA_BUCKET_NAME"
fi

# Get most recent Cypher scripts
echo "Fetching most recent Cypher scripts"
aws s3 cp --recursive s3://$DATA_BUCKET_NAME/$NEO4J_CYPHER_PATH_S3/ $NEO4J_CYPHER_PATH

# Download data to /var/lib/neo4j/import
echo "Downloading CSV data for release $RELEASE"
aws s3 cp --recursive s3://$DATA_BUCKET_NAME/data/$RELEASE/csv/ $NEO4J_IMPORT_PATH/

# Update Cypher load query for correct release
mkdir -p $NEO4J_CYPHER_PATH/tmp/$RELEASE/
cat /var/lib/neo4j/cypher/load.cyp | sed "s/RELEASE/$RELEASE/g" > $NEO4J_CYPHER_PATH/tmp/$RELEASE/load.$RELEASE.cyp

printf "Updated script for release $RELEASE:\n$(cat $NEO4J_CYPHER_PATH/tmp/$RELEASE/load.$RELEASE.cyp)"

# Run Cypher load query
echo "Loading data for release $RELEASE into Neo4j..."
cat $NEO4J_CYPHER_PATH/tmp/$RELEASE/load.$RELEASE.cyp | \
    /usr/bin/cypher-shell \
        --username $NEO4J_USERNAME \
        --password $NEO4J_PASSWORD \
        --format verbose

# TODO: Conditional queries for alignments, KIR (requires running separate Cypher scripts)
# if $ALIGN; then \
    # load alignments

echo "Done"
exit 0