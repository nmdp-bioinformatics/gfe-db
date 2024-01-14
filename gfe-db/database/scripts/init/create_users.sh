#!/bin/bash -xeu

# Comma-separated list of users to create
USERS=$1 # username:password,...

if [ -z "${USERS}" ]; then
    echo "No users specified, exiting"
    exit 1
fi

if [ "${USE_PRIVATE_SUBNET}" = true ]; then
  NEO4J_URI=bolt://127.0.0.1:7687
  ENCRYPTION=false
else
    NEO4J_URI=neo4j+s://${SUBDOMAIN}.${HOST_DOMAIN}:7687
    ENCRYPTION=default
fi

echo "NEO4J_URI=${NEO4J_URI}"
echo "ENCRYPTION=${ENCRYPTION}"

RETURN_CODES=0
for user in $(echo $USERS | sed "s/,/ /g")
do
    username=$(echo $user | cut -d':' -f1)
    password=$(echo $user | cut -d':' -f2)
    echo "Creating user ${username}..."
    cat $NEO4J_HOME/cypher/create_user.cyp | \
    cypher-shell \
        -u neo4j \
        -p ${NEO4J_PASSWORD} \
        -a ${NEO4J_URI} \
        --encryption ${ENCRYPTION} \
        -P "username => \"${username}\"" \
        -P "password => \"${password}\""
    RETURN_CODES=$((RETURN_CODES + $?))
done

if [ $RETURN_CODES -ne 0 ]; then
    echo "Failed to create users"
    exit 1
else
    echo "Success"
    exit 0
fi