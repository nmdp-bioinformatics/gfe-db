#!/bin/bash

export ROOT=$(dirname $(dirname "$0"))
export SRC_DIR=$ROOT/src
export CYPHER_DIR=$ROOT/cypher
export LOAD_SCRIPT=load.cyp

# TODO: Retrieve NEO4J_HOST from SSM Parameter store
# Check for environment variables
if [[ "${NEO4J_HOST_SSM_PARAM}" ]]; then
	export NEO4J_HOST="$(aws ssm get-parameter \
		--name "$NEO4J_HOST_SSM_PARAM" \
		--region $REGION \
	| jq -r '.Parameter.Value')"
	echo "Found NEO4J_HOST_SSM_PARAM: $NEO4J_HOST_SSM_PARAM"
else
	echo "NEO4J_HOST_SSM_PARAM not set. Looking for NEO4J_HOST..."
fi

if [[ -z "${NEO4J_HOST}" ]]; then
	echo "NEO4J_HOST not set"
	exit 1
elif [[ -z "${GFE_BUCKET}" ]]; then
	echo "GFE_BUCKET not set"
	exit 1
elif [[ -z "${RELEASES}" ]]; then
	echo "RELEASES not set. Please specify the release versions to load."
	exit 1
elif [[ -z "${NEO4J_USERNAME}" ]]; then
	echo "NEO4J_USERNAME not set"
	exit 1
elif [[ -z "${NEO4J_PASSWORD}" ]]; then
	echo "NEO4J_PASSWORD not set"
	exit 1
fi

echo -e "Found environment variables:"
echo -e "GFE_BUCKET: $GFE_BUCKET\n\
RELEASES: $RELEASES\n\
NEO4J_HOST_SSM_PARAM: $NEO4J_HOST_SSM_PARAM\n\
NEO4J_HOST: $NEO4J_HOST\n\
NEO4J_USERNAME: $NEO4J_USERNAME\n\
NEO4J_PASSWORD: $NEO4J_PASSWORD"

# Load Neo4j through the HTTP API
echo
echo "Loading Neo4j database..."
python3 $SRC_DIR/load.py
