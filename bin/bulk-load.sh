#!/bin/sh

# set this to where the csv and graph will live
DATA_HOME=`pwd`

# path to use from laptop
export NEO4J_HOME=/var/lib/neo4j

#IMPORTCMD=$NEO4J_HOME/bin/neo4j-admin import
IMPORTCMD=$NEO4J_HOME/bin/neo4j-import
NEO4JCMD=$NEO4J_HOME/bin/neo4j
GRAPHPATH=/var/lib/neo4j/graph/databases/imputewmda.db

$NEO4JCMD stop

rm -rf $GRAPHPATH/*


# # import graph
set -x
$IMPORTCMD --into $GRAPHPATH \
	--id-type INTEGER \
	--nodes /opt/allele_Nodes.csv \
	--nodes /opt/sequence_Nodes.csv \
    --relationships /opt/group_edges.csv \
    --relationships /opt/gfe_edges.csv \
    --relationships /opt/sequence_edges.csv  


# start neo4j
$NEO4JCMD start
