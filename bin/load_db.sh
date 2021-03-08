#!/bin/bash

CYPHER=load_duplicate_rels.cyp

echo "Load script: $CYPHER"
# TO DO: format IMGTHLA version string
for RELEASE in $(echo $RELEASES | sed "s/,/ /g")
do
    echo "Loading release version $RELEASE..."
    cat neo4j/$CYPHER | \
        sed "s/RELEASE/$RELEASE/g" | \
        docker exec --interactive gfe cypher-shell -u neo4j -p gfedb
done

echo "Done"