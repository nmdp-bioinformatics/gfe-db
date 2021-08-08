#!/bin/bash

export ROOT=$(dirname $(dirname "$0"))
export SRC_DIR=$ROOT/src
export CYPHER_DIR=$ROOT/cypher
export LOAD_SCRIPT=load.cyp

export AWS_ACCESS_KEY_ID="ASIAXXVOWRIZGZGVD6W4"
export AWS_SECRET_ACCESS_KEY="p7+h689PJLzte/ofhMrR5kN/KdnOvI+Z3J24gsa8"
export AWS_SESSION_TOKEN="IQoJb3JpZ2luX2VjEGoaCXVzLWVhc3QtMSJHMEUCIQCJoY/93NpJjDmXoPH0+Rg8IWQ9gCboePHhix7IdNzxVAIgPFmhsf1eoPa5DbsbS7LSZEcL87V8vFdYgeoTWT9rMhAqgQMIg///////////ARACGgw1MzE4Njg1ODQ0OTgiDCy1XI/U1Q07HxjPvCrVApGunuG+tTB70D1lTTCOUUSeHA5EmQqWb1G5FdpeW+pzzMmaFCx2lYNk9tFXYWhhDPWuIBwp7jtAcVFvGi2R4yRFAHhTB+7fCS8Zu+CptMybmcHgviq+YQ1GedOjk9F+SwzNmGVcQHp9PD2MxElxW5podtvbaZMxeYUQl9VIwPQiDai1DQwANJHDpL+TzHDVK5MB9h7YQ14rapg9dgnQqaN6NbUpw45IcNq/iswp/GEmSNS7Lv3Yw0QDN3GT6THkInDVPtlGuYz/pyQB1t+H2hz9jo5zqCfghaQa1LG+SH3EKxf0IRIeXgGbsfMs16SmLiAjNsSYQCWt73+QgGX+LQu/xvB2UeIwTLn7mMvFi69q2ExwupoffrpWhvQMA4zFqT47s0qsfUFIOwqOZJvegOPLlYJClb6hajLLwkoKQNV/STAmhTMrsCa+mV4pvfhPrgMEl9ZOMJDovIgGOqYBdYUcUTLIij2RVOv6yHPEXZw6L2bsqtTzPyw41GYjZR1BUfqTIas73tHvNeq1L5/fUvl3O5/EkIGxtgvcrP5TurkAT/cj+C6HQGYzvER2D1if6XR3Wdn3m3FtlmLafcQ+bLERzKrgsS4jAyvwBOn5GOpVTWCI6NDTJjdPIrg85TY+YaA0mDR4Rq78KkQZGi2N06mixtIu8I/bfz/PUdB88LqDtFf3bw=="

NEO4J_HOST=44.192.54.30
NEO4J_USERNAME=neo4j
NEO4J_PASSWORD=gfedb

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
echo "Loading Neo4j database..."
python3 $SRC_DIR/load_gfedb.py
