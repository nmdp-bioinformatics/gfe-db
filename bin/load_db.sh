#!/bin/bash

for RELEASE in $(echo $RELEASES | sed "s/,/ /g")
do
    cat neo4j/load_new.cyp | \
        sed "s/RELEASE/$RELEASE/g" | \
        docker exec --interactive gfe cypher-shell -u neo4j -p gfedb
done