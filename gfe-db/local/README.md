# Local Deployment for gfe-db
Instructions for using Docker to run gfe-db locally.

## Features
- Set password
- Persists data locally so the container can be stopped and started

## Usage
Run a separate container before running the Neo4j container to prepare the data from the backup.
```bash
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
```

Basic run command.
```bash
docker run \
    --restart always \
    --publish=7474:7474 --publish=7687:7687 \
    --env NEO4J_AUTH=neo4j/gfedb2023 \
    --volume=$(pwd)/neo4j/data:/data \
    --volume=$(pwd)/neo4j/logs:/logs \
    neo4j:5.15
```

docker run \
    --restart always \
    --publish=7474:7474 --publish=7687:7687 \
    --env NEO4J_AUTH=neo4j/gfedb2023 \
    --volume=$(pwd)/neo4j/data:/data \
    --volume=$HOME/neo4j/logs:/logs \
    neo4j:5.15

## Development Steps
- 