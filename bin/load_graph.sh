#!/usr/bin/env bash

# Set JAVA_HOME

# Script variables
BIN=$(dirname "$0")

if [ "$JAVA_HOME" ]; then
	echo "Java Home is set at " ${JAVA_HOME}
else
	echo "Exiting. JAVA_HOME env variable is not set."
	exit
fi

# Set NEO4J_HOME 
if [ "$NEO4J_HOME" ]; then
	echo "Neo4j is at " ${NEO4J_HOME}
else
	echo "Exiting. NEO4J_HOME env variable is not set."
	exit
fi

NEO4J_IMPORT_CMD=$NEO4J_HOME/bin/neo4j-import
NEO4J_CMD=$NEO4J_HOME/bin/neo4j

# DATA_HOME will have the csv output and graph database and config dir
#DATA_HOME=${PWD}/output
DATA_HOME=$1


NEO4J_ADMIN_CMD=${NEO4J_HOME}/bin/neo4j-admin
NEO4J_DEFAULT_PASSWORD=gfedb

NEO4J_ACTIVE_DATABASE=gfedb.db
GRAPH_PATH=${DATA_HOME}/databases/${NEO4J_ACTIVE_DATABASE}

#IMPORT_PATH=/opt
IMPORT_PATH=${DATA_HOME}

echo $IMPORT_PATH

${NEO4J_CMD} stop

# Remove the previous database
rm -rf ${GRAPH_PATH}/*

# Copy the config to the conf directory
#cp neo4j/conf/neo4j.conf.template ${NEO4J_CONF_DIR}/neo4j.conf
#cp ${BIN}/neo4j.conf.template ${NEO4J_CONF_DIR}/neo4j.conf
cp /opt/neo4j.conf.template /var/lib/neo4j/conf/neo4j.conf

# Create conf directory
NEO4J_CONF_DIR=${DATA_HOME}/conf
mkdir -p ${NEO4J_CONF_DIR}

# Copy the config to the conf directory
#cp neo4j/conf/neo4j.conf.template ${NEO4J_CONF_DIR}/neo4j.conf
cp ${BIN}/neo4j.conf.template ${NEO4J_CONF_DIR}/neo4j.conf

# Update neo4j.conf file
echo dbms.directories.data=${DATA_HOME} >> ${NEO4J_CONF_DIR}/neo4j.conf
echo dbms.active_database=${NEO4J_ACTIVE_DATABASE} >> ${NEO4J_CONF_DIR}/neo4j.conf

# import CVS into the graph database
${NEO4J_IMPORT_CMD} --into ${GRAPH_PATH} \
	--id-type INTEGER \
	--nodes ${IMPORT_PATH}/allele_nodes.csv \
	--nodes ${IMPORT_PATH}/sequence_nodes.csv \
	--nodes ${IMPORT_PATH}/cds_nodes.csv \
    --relationships ${IMPORT_PATH}/gfe_edges.csv \
    --relationships ${IMPORT_PATH}/seq_edges.csv  \
    --relationships ${IMPORT_PATH}/group_edges.csv \
    --relationships ${IMPORT_PATH}/cds_edges.csv \


# Set default password to $NEO4J_DEFAULT_PASSWORD
rm -f output/dbms/auth
NEO4J_CONF=${NEO4J_CONF_DIR} ${NEO4J_ADMIN_CMD} set-initial-password ${NEO4J_DEFAULT_PASSWORD}

# start neo4j
NEO4J_CONF=${NEO4J_CONF_DIR} ${NEO4J_CMD} console


