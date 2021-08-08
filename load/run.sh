#!/bin/bash

export ROOT=$(dirname $(dirname "$0"))
export SRC_DIR=$ROOT/load/src
export CYPHER_DIR=$ROOT/load/cypher
export LOAD_SCRIPT=load.cyp

# Check if RELEASES is set
if [ -z ${RELEASES+x} ]; then 
    echo "RELEASES is not set. Please specify the release versions to load."; 
else 
    echo "Loading IMGT/HLA release versions: $RELEASES";
fi

# Load Neo4j through the HTTP API
python3 $SRC_DIR/load_gfedb.py
