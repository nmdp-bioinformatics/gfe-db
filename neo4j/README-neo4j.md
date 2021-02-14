# Running Neo4j 4.2 in Docker

## Build the image
```
docker build --tag gfe-db .
```

## Run the container
```
docker run -p 7474:7474 -p 7473:7473 -p 7687:7687 gfe-db
```