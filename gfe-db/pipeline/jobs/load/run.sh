#!/bin/bash

export ROOT=$(dirname $(dirname "$0"))
export SRC_DIR=$ROOT/src
export CYPHER_DIR=$ROOT/cypher
export LOAD_SCRIPT=load.cyp

# Check for environment variables
if [[ -z "${GFE_BUCKET}" ]]; then
	echo "GFE_BUCKET not set"
	exit 1
elif [[ -z "${RELEASES}" ]]; then
	echo "RELEASES not set. Please specify the release versions to load."
	exit 1
elif [[ -z "${NEO4J_HOST}" ]]; then
	echo "NEO4J_HOST not set"
	exit 1
elif [[ -z "${NEO4J_USERNAME}" ]]; then
	echo "NEO4J_USERNAME not set"
	exit 1
elif [[ -z "${NEO4J_PASSWORD}" ]]; then
	echo "NEO4J_PASSWORD not set"
	exit 1
else
	echo "Found environment variables:"
	echo -e "GFE_BUCKET: $GFE_BUCKET\nRELEASES: $RELEASES\nNEO4J_HOST: $NEO4J_HOST\nNEO4J_USERNAME: $NEO4J_USERNAME\nNEO4J_PASSWORD: $NEO4J_PASSWORD"
fi

# Load Neo4j through the HTTP API
echo
echo "Loading Neo4j database..."
python3 $SRC_DIR/load_gfedb.py
