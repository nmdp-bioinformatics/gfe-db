#!/bin/bash

# CYPHER_PATH=neo4j/cypher
# # SCRIPT=load.cyp
# SCRIPT=test.cyp

PROTOCOL=http
HOST=35.173.36.115

# NEO4J_USER="neo4j"
# NEO4J_PASSWORD="gfedb"
PORT=7474
ENDPOINT=db/neo4j/tx/commit
URL=$PROTOCOL://$HOST:$PORT/$ENDPOINT

echo "Load script: $CYPHER_PATH/$SCRIPT"
echo "URL: $URL"

# Check if RELEASES is set
if [ -z ${RELEASES+x} ]; then 
    echo "RELEASES is not set. Please specify the release versions to load."; 
else 
    echo "Loading IMGT/HLA release versions: $RELEASES";
fi

echo "Creating constraints and indexes..."
cat $CYPHER_PATH/create_index.cyp | \
    # docker exec -i gfe cypher-shell -u neo4j -p gfedb

# TO DO: format IMGTHLA version string
for RELEASE in $(echo $RELEASES | sed "s/,/ /g"); do
    echo "Loading release version $RELEASE..."
    start_release=$SECONDS
    # cat $CYPHER_PATH/$SCRIPT | \
    #     sed "s/RELEASE/$RELEASE/g" | \
    #     docker exec -i gfe cypher-shell -u neo4j -p gfedb

    # Load via Neo4j HTTP API
    statement="\"$(cat $CYPHER_PATH/$SCRIPT)\""
    payload='{"statements": [{"statement": '$statement', "params": {}}]}'
    # echo $statement
    echo "Username/password: $NEO4J_USER:$NEO4J_PASSWORD"
    echo $(echo -n $NEO4J_USER:$NEO4J_PASSWORD | base64)
    echo $payload | \
        curl \
            --data-binary @- \
            -H "content-type:application/json" \
            -H "accept:application/json" \
            -H "authorization: Basic bmVvNGo6Z2ZlZGI=" \
            -X POST $URL

    duration=$(( SECONDS - start_release ))
    echo "Loaded in $duration seconds"
    # echo -e "\n"
done
echo "Done"

