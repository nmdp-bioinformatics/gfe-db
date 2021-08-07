#!/bin/bash

export ROOT=$(dirname $(dirname "$0"))
export BIN_DIR=$ROOT/scripts
export SRC_DIR=$ROOT/src
export DATA_DIR=$ROOT/data
export LOGS_DIR=$ROOT/logs
export CYPHER_PATH=neo4j/cypher
export SCRIPT=load.cyp

export GFE_BUCKET=gfe-db-4498
export RELEASES="3440"
export ALIGN=True
export KIR=False
export MEM_PROFILE=True
# export CYPHER_PATH=neo4j/cypher
# export SCRIPT=load.cyp