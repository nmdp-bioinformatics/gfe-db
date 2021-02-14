# Running Neo4j 4.2 in Docker

## Build the HLA data
```
bash bin/buld_gfedb.py
```

## Build Neo4j Docker image
```
docker build --tag gfe-db .
```

## Run the container
```
docker run -p 7474:7474 -p 7473:7473 -p 7687:7687 gfe-db
```

## Populate the Graph
<!-- TO DO: figure out shell command to load graph -->