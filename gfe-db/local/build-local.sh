#!/bin/bash -eu

# TODO create a version tag from the data path for Docker, push both latest and versioned tags
DATA_S3_PATH=$1
TEMP_DIR=/tmp/$APP_NAME

# TODO validate the path using the same logic as the Makefile

# Create volumes to mount
mkdir -p $TEMP_DIR neo4j/{data,logs,backups/{system,neo4j}}

# Download most recent backup
aws s3 cp $DATA_S3_PATH $TEMP_DIR/gfedb.zip
unzip -o $TEMP_DIR/gfedb.zip -d neo4j

# Load backup into volumes
docker run --interactive --tty --rm \
    --platform=linux/amd64 \
    --volume=$(pwd)/neo4j/data:/data \
    --volume=$(pwd)/neo4j/backups:/backups \
    neo4j/neo4j-admin:5.15.0 \
    neo4j-admin database load --from-path=/backups/system --overwrite-destination=true --verbose system

docker run --interactive --tty --rm \
    --platform=linux/amd64 \
    --volume=$(pwd)/neo4j/data:/data \
    --volume=$(pwd)/neo4j/backups:/backups \
    neo4j/neo4j-admin:5.15.0 \
    neo4j-admin database load --from-path=/backups/neo4j --overwrite-destination=true --verbose neo4j

# Build, tag and push image
docker build -t gfe-db .
docker tag gfe-db:latest "$DOCKER_USERNAME/gfe-db:latest"
docker login -u $DOCKER_USERNAME -p $DOCKER_PASSWORD
docker push "$DOCKER_USERNAME/gfe-db:latest"

# Clean up
rm -rf neo4j/{data,backups} /tmp/$APP_NAME
