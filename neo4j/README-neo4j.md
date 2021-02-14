# Running Neo4j 4.2 in Docker

## Build the image
```
docker build --file /neo4j/Dockerfile --tag gfedb-neo4j-4
```

## Run the container
```
docker run --interactive --tty --name custom-container-1 -p7687:7687 -p7474:7474 -p7473:7473 --env NEO4J_AUTH=neo4j/password --env NEO4J_ACCEPT_LICENSE_AGREEMENT=yes neo4j:4.2.3-enterprise-custom-container-1
```
```
docker run --name neo4j-gfe-db -d -p 7474:7474 -p 7473:7473 -p 7687:7687 gfedb-neo4j-4
```