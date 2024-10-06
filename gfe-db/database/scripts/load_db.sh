#!/bin/bash -x

source /home/ec2-user/env.sh

# Get APP_NAME, AWS_REGION, STAGE setup on db install
if [ -z $EC2_USER_HOME ]; then
    echo "ERROR: EC2_USER_HOME not set"
    exit 1
fi
if [ -z $APP_NAME ]; then
    echo "ERROR: APP_NAME not set"
    exit 1
fi

if [ -z $STAGE ]; then
    echo "ERROR: STAGE not set"
    exit 1
fi

if [ -z $SERVICE_NAME ]; then
    echo "ERROR: SERVICE_NAME not set"
    exit 1
fi

if [[ -z $NEO4J_HOME ]]; then
    echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Neo4j not found"
    exit 1
else
    echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Found Neo4j in $NEO4J_HOME"
fi

# Check for release argument
RELEASE=$1

# Set paths
NEO4J_CYPHER_PATH=$NEO4J_HOME/cypher
NEO4J_IMPORT_PATH=$NEO4J_HOME/import
S3_NEO4J_CYPHER_PATH=config/$SERVICE_NAME/neo4j/cypher

# TODO Get from state payload
S3_CSV_PATH=data/$RELEASE/csv

if [[ -z $AWS_REGION ]]; then
    export AWS_REGION=$(curl --silent http://169.254.169.254/latest/dynamic/instance-identity/document | jq -r '.region')
fi

if [[ -z $RELEASE ]]; then
    echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Release version not found"
    exit 1
else
    echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Starting load process for $RELEASE"
fi

# Get Neo4j Credentials
NEO4J_CREDENTIALS=$(aws secretsmanager get-secret-value \
    --region $AWS_REGION \
    --secret-id /${APP_NAME}/${STAGE}/${AWS_REGION}/Neo4jCredentials | jq -r '.SecretString')
NEO4J_USERNAME=$(echo $NEO4J_CREDENTIALS | jq -r '.NEO4J_USERNAME')
NEO4J_PASSWORD=$(echo $NEO4J_CREDENTIALS | jq -r '.NEO4J_PASSWORD')

# Get data bucket name
DATA_BUCKET_NAME=$(aws ssm get-parameters \
    --region $(curl -s http://169.254.169.254/latest/meta-data/placement/availability-zone | sed 's/\(.*\)[a-z]/\1/') \
    --names "/${APP_NAME}/${STAGE}/${AWS_REGION}/DataBucketName" \
    | jq -r '.Parameters[0].Value')

if [[ -z $DATA_BUCKET_NAME ]]; then
    echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - S3 bucket not found."
    exit 1
else
    echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Found S3 bucket: $DATA_BUCKET_NAME"
fi

# Get most recent Cypher scripts
echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Fetching most recent Cypher scripts"
aws s3 cp --recursive s3://$DATA_BUCKET_NAME/$S3_NEO4J_CYPHER_PATH/ $NEO4J_CYPHER_PATH --quiet
# check error status of aws s3 cp and abort if not zero
[ $? -eq 0 ] || exit 1
# TODO validate file was downloaded, abort if not, so that a failure signal can be sent to Step Functions

# Download data to NEO4J_HOME/import
echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Downloading CSV data for release $RELEASE"
aws s3 cp --recursive s3://$DATA_BUCKET_NAME/$S3_CSV_PATH/ $NEO4J_IMPORT_PATH/ --quiet

# Update Cypher load query for correct release
# TODO Change load.cyp to load.cyp.template
mkdir -p $NEO4J_CYPHER_PATH/tmp/$RELEASE/

# TODO Use Cypher params for RELEASE instead of sed
cat $NEO4J_CYPHER_PATH/load.cyp | sed "s/RELEASE/$RELEASE/g" > $NEO4J_CYPHER_PATH/tmp/$RELEASE/load.$RELEASE.cyp

# check error status of sed and abort if not zero
[ -f $NEO4J_CYPHER_PATH/tmp/$RELEASE/load.$RELEASE.cyp ] || exit 1

echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Executing query"
echo "****** Begin Cypher ******"
printf "$(cat $NEO4J_CYPHER_PATH/tmp/$RELEASE/load.$RELEASE.cyp)\n"
echo "****** End Cypher ******"

# Run Cypher load query
echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Loading data for release $RELEASE into Neo4j..."

if [[ "$USE_PRIVATE_SUBNET" = true ]]; then

    # # With SSL/TLS policy disabled for private instance
    cat $NEO4J_CYPHER_PATH/tmp/$RELEASE/load.$RELEASE.cyp | \
        $NEO4J_HOME/bin/cypher-shell \
            --address bolt://127.0.0.1:7687 \
            --encryption false \
            --username $NEO4J_USERNAME \
            --password $NEO4J_PASSWORD \
            --format verbose
    LOAD_EXIT_STATUS=$?

else

    # With SSL/TLS policy enabled
    cat $NEO4J_CYPHER_PATH/tmp/$RELEASE/load.$RELEASE.cyp | \
        $NEO4J_HOME/bin/cypher-shell \
            --address neo4j://$SUBDOMAIN.$HOST_DOMAIN:7687 \
            --encryption true \
            --username $NEO4J_USERNAME \
            --password $NEO4J_PASSWORD \
            --format verbose
    LOAD_EXIT_STATUS=$?

fi

if [[ $LOAD_EXIT_STATUS -eq 0 ]]; then
    echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Load complete"
else
    echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Load failed"
fi

# TODO: Conditional queries for alignments, KIR (requires running separate Cypher scripts)
# if $ALIGN; then \
    # load alignments

# TODO: if $? == 0 for all queries, send TaskSuccess to StepFunctions API
echo "$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Cleaning up"
rm -r $NEO4J_IMPORT_PATH/*

exit $LOAD_EXIT_STATUS
