# Running Neo4j 4.2 in Docker

## Build the HLA data
```
bash bin/build.sh
```

## Build Neo4j Docker image
```
docker build --tag gfe-db .
```

## Run the container
```
docker run --name db -v "$(pwd)"/data/csv/:/var/lib/neo4j/import \
    -p 7474:7474 -p 7473:7473 -p 7687:7687 gfe-db
```

## Populate the Graph
Once the Neo4j server is running, the data can be loaded using the cypher script:
```
cat neo4j/update_hla.cyp | docker exec --interactive db cypher-shell -u neo4j -p gfedb
```