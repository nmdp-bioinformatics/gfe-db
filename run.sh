#!/bin/bash

neo4j/gather_neo4j_plugins.sh
docker build --tag gfe-db . && \
docker run -p 7474:7474 -p 7473:7473 -p 7687:7687 gfe-db