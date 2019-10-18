#!/usr/bin/env bash

#IMPORT_PATH will have the csv output
IMPORT_PATH=${1}

NEO4J_DEFAULT_PASSWORD=gfedb

# The default graph database
mkdir -p /var/lib/neo4j/gfedb/
GRAPH_PATH=/var/lib/neo4j/gfedb/graph.db

echo Importing CSV files from $IMPORT_PATH

neo4j stop

# import CVS into the graph database
neo4j-import --into ${GRAPH_PATH} \
	--id-type INTEGER \
	--nodes ${IMPORT_PATH}/allele_nodes.csv \
	--nodes ${IMPORT_PATH}/sequence_nodes.csv \
	--nodes ${IMPORT_PATH}/cds_nodes.csv \
    --relationships ${IMPORT_PATH}/gfe_edges.csv \
    --relationships ${IMPORT_PATH}/seq_edges.csv  \
    --relationships ${IMPORT_PATH}/group_edges.csv \
    --relationships ${IMPORT_PATH}/cds_edges.csv 


# Set default password to $NEO4J_DEFAULT_PASSWORD
neo4j-admin set-initial-password ${NEO4J_DEFAULT_PASSWORD}

