FROM neo4j:4.2

RUN apt-get update \
    && apt-get install -y curl openssl apt-utils zip unzip

ENV NEO4J_AUTH=neo4j/gfedb \
    NEO4J_ACCEPT_LICENSE_AGREEMENT=yes \
    NEO4J_dbms_memory_heap_initial__size=2G \
    NEO4J_dbms_memory_heap_max__size=2G \
    NEO4J_dbms_memory_pagecache_size=1G \
    NEO4J_dbms_security_procedures_unrestricted=apoc.*,gds.* \
    NEO4J_dbms_security_allow__csv__import__from__file__urls=true \
    # NEO4J_apoc_import_file_enabled=true \
    # NEO4J_apoc_import_file_use__neo4j__config=true \
    # NEO4J_apoc_export_file_enabled=true 
    GDS_LIB_VERSION=1.4.1 \ 
    APOC_LIB_VERSION=4.2.0.1 \
    GITHUB_GDS_URI=https://github.com/neo4j/graph-data-science/releases/download \
    GITHUB_APOC_URI=https://github.com/neo4j-contrib/neo4j-apoc-procedures/releases/download \
    NEO4J_GDS_URI=${GITHUB_GDS_URI}/${GDS_LIB_VERSION}/neo4j-graph-data-science-${GDS_LIB_VERSION}-standalone.jar \
    NEO4J_APOC_URI=${GITHUB_APOC_URI}/${APOC_LIB_VERSION}/apoc-${APOC_LIB_VERSION}-all.jar

# Download Neo4j APOC libraries
RUN sh -c 'cd /var/lib/neo4j/plugins && \
    echo "Downloading APOC libraries..." \
    curl -C- --progress-bar \
        --location ${NEO4J_APOC_URI} \
        --output $NEO4J_DIR/plugins/apoc-${APOC_LIB_VERSION}-all.jar'

# Download Neo4j GDS libraries
RUN sh -c 'cd /var/lib/neo4j/plugins && \
    echo "Downloading Neo4j Graph Data Science libraries..." \
    curl -C- --progress-bar \
        --location ${NEO4J_GDS_URI} \
        --output $NEO4J_DIR/plugins/neo4j-graph-data-science-${GDS_LIB_VERSION}-standalone.jar'

# Copy GFE data as CSV files
# COPY data/csv/ /var/lib/neo4j/import/

# Mount the data directory directly into Neo4j import directory
VOLUME /data/csv/ /var/lib/neo4j/import/

EXPOSE 7474 7473 7687

CMD ["neo4j"]
