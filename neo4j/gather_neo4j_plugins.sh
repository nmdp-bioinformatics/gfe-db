#!/bin/bash

####################
# Download for neo4j into plugins/ directory following libraries
# - neo4j-apoc-procedures
# - graph-data-science
####################

NEO4J_DIR=neo4j
NEO4J_VERSION=4.2.2
GDS_LIB_VERSION=1.4.1
APOC_LIB_VERSION=4.2.0.1
GITHUB_GDS_URI=https://github.com/neo4j/graph-data-science/releases/download
GITHUB_APOC_URI=https://github.com/neo4j-contrib/neo4j-apoc-procedures/releases/download
NEO4J_GDS_URI=${GITHUB_GDS_URI}/${GDS_LIB_VERSION}/neo4j-graph-data-science-${GDS_LIB_VERSION}-standalone.jar
NEO4J_APOC_URI=${GITHUB_APOC_URI}/${APOC_LIB_VERSION}/apoc-${APOC_LIB_VERSION}-all.jar

mkdir -p $NEO4J_DIR/plugins

echo "Downloading Neo4j Graph Data Science libraries..."
curl -C- --progress-bar \
    --location ${NEO4J_GDS_URI} \
    --output $NEO4J_DIR/plugins/neo4j-graph-data-science-${GDS_LIB_VERSION}-standalone.jar
echo "Downloading APOC libraries..."
curl -C- --progress-bar \
    --location ${NEO4J_APOC_URI} \
    --output $NEO4J_DIR/plugins/apoc-${APOC_LIB_VERSION}-all.jar
echo 'Done.'
