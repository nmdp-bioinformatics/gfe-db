#!/bin/bash

export ROOT=$(pwd)
export LOAD_DIR=$ROOT/load
export NEO4J_DIR=$ROOT/neo4j
export SCRIPT=load.cyp

# Check if RELEASES is set
if [ -z ${RELEASES+x} ]; then 
    echo "RELEASES is not set. Please specify the release versions to load."; 
else 
    echo "Loading IMGT/HLA release versions: $RELEASES";
fi

# Load Neo4j thruogh the HTTP API
python3 "$LOAD_DIR"/src/load_gfedb.py
