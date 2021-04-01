#!/bin/bash

CYPHER_PATH=neo4j
SCRIPT=load.cyp

echo "Load script: $CYPHER_PATH/$SCRIPT"

# echo "Creating constraints and indexes..."
# cat $CYPHER_PATH/create_index.cyp | \
#     docker exec -i gfe cypher-shell -u neo4j -p gfedb

# TO DO: format IMGTHLA version string
for RELEASE in $(echo $RELEASES | sed "s/,/ /g")
do
    echo "Loading release version $RELEASE..."
    start_release=$SECONDS
    cat $CYPHER_PATH/$SCRIPT | \
        sed "s/RELEASE/$RELEASE/g" | \
        docker exec -i gfe cypher-shell -u neo4j -p gfedb
    duration=$(( SECONDS - start_release ))
    echo "$duration seconds"
done
echo "Done"