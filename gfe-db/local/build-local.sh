#!/bin/bash -eu

DATA_S3_PATH=$1
TEMP_DIR=/tmp/$APP_NAME

# Check that NEO4J_DATABASE_NAME is set
if [ -z ${NEO4J_DATABASE_NAME+x} ]; then
    echo "NEO4J_DATABASE_NAME is not set" >&2
    exit 1
else
    echo "NEO4J_DATABASE_NAME is set to '$NEO4J_DATABASE_NAME'"
fi

format_s3_date() {
    local s3_path="$1"

    # Regular expression to extract date and time parts
    if [[ $s3_path =~ ([0-9]{4})/([0-9]{2})/([0-9]{2})/([0-9]{2})/([0-9]{2}) ]]; then
        local year=${BASH_REMATCH[1]}
        local month=${BASH_REMATCH[2]}
        local day=${BASH_REMATCH[3]}
        local hour=${BASH_REMATCH[4]}
        local minute=${BASH_REMATCH[5]}

        # Combine the parts into the desired format
        echo "${year}${month}${day}_${hour}${minute}"
    else
        echo "Could not parse date from S3 path: $s3_path" >&2
        exit 1
    fi
}

echo "This version is $(format_s3_date $DATA_S3_PATH)"

# Create volumes to mount
mkdir -p $TEMP_DIR neo4j/{data,logs,backups/{system,$NEO4J_DATABASE_NAME}}

# Download most recent backup
aws s3 cp $DATA_S3_PATH $TEMP_DIR/$NEO4J_DATABASE_NAME.zip
unzip -o $TEMP_DIR/$NEO4J_DATABASE_NAME.zip -d neo4j

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
    neo4j-admin database load --from-path=/backups/$NEO4J_DATABASE_NAME --overwrite-destination=true --verbose $NEO4J_DATABASE_NAME

# Build, tag and push image
version_tag=$(format_s3_date $DATA_S3_PATH)
docker build --build-arg NEO4J_DATABASE_NAME=$NEO4J_DATABASE_NAME -t gfe-db .
docker tag gfe-db:latest $DOCKER_USERNAME/gfe-db:latest
docker tag gfe-db:latest $DOCKER_USERNAME/gfe-db:$version_tag
echo $DOCKER_PASSWORD | docker login -u $DOCKER_USERNAME --password-stdin
docker push $DOCKER_USERNAME/gfe-db:latest
docker push $DOCKER_USERNAME/gfe-db:$version_tag

# Clean up
rm -rf neo4j/{data,backups} /tmp/$APP_NAME
