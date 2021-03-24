#!/bin/bash

CYPHER_PATH=neo4j/load_test.cyp

echo "Load script: $CYPHER_PATH"
# TO DO: format IMGTHLA version string
for RELEASE in $(echo $RELEASES | sed "s/,/ /g")
do
    echo "Loading release version $RELEASE..."
    cat $CYPHER_PATH | \
        sed "s/RELEASE/$RELEASE/g" | \
        docker exec -i gfe cypher-shell -u neo4j -p gfedb
done

echo "Done"