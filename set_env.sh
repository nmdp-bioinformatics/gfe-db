#!/bin/bash

export ROOT=$(pwd) # $(dirname $(dirname "$0"))
export BIN_DIR=$ROOT/scripts
export SRC_DIR=$ROOT/src
export DATA_DIR=$ROOT/data
export LOGS_DIR=$ROOT/logs
export NEO4J_DIR=$ROOT/neo4j
export SCRIPT=load.cyp

export GFE_BUCKET=gfe-db-4498
export RELEASES="3440"
export ALIGN=True
export KIR=False
export MEM_PROFILE=True

export NEO4J_HOST=44.192.54.30
export NEO4J_USERNAME=neo4j
export NEO4J_PASSWORD=gfedb