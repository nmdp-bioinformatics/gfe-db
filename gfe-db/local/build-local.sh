#!/bin/bash -eu

from_path=s3://nmdpf-gfe-db-810526023897-us-east-1/backups/neo4j/2024/01/07/20/18/gfedb.zip # $1

# Create volumes to mount
mkdir -p /tmp/restore neo4j/{data,logs,backups/{system,neo4j}}

# Download most recent backup
aws s3 cp $from_path /tmp/restore/gfedb.zip
unzip -o /tmp/restore/gfedb.zip -d neo4j

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

# Build image
docker build -t gfe-db .