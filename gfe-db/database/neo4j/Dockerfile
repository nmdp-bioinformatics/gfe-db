FROM neo4j:latest

RUN apt-get update
    # \ && apt-get install -y curl openssl apt-utils zip unzip

ENV NEO4J_AUTH=neo4j/gfedb
ENV NEO4J_ACCEPT_LICENSE_AGREEMENT=yes
ENV NEO4J_dbms_memory_heap_initial__size=8G
ENV NEO4J_dbms_memory_heap_max__size=8G
ENV NEO4J_dbms_memory_pagecache_size=4G
# ENV NEO4J_dbms_memory_transaction_global__max__size=1G
ENV NEO4J_dbms_security_auth__enabled=true
ENV NEO4J_dbms_security_procedures_unrestricted=apoc.*,gds.*
ENV NEO4J_dbms_security_allow__csv__import__from__file__urls=true
ENV NEO4J_apoc_import_file_enabled=true
ENV NEO4J_apoc_import_file_use__neo4j__config=true
ENV NEO4J_apoc_export_file_enabled=true 
ENV NEO4J_dbms_default__listen__address=0.0.0.0
# ENV NEO4J_dbms_connector_https_enabled=true
# ENV NEO4J_dbms_logs_user_stdout__enabled=true
# ENV NEO4JLABS_PLUGINS='["apoc","graph-data-science"]'
# ENV GDS_LIB_VERSION=1.4.1
# ENV APOC_LIB_VERSION=4.2.0.2
# ENV GITHUB_GDS_URI=https://github.com/neo4j/graph-data-science/releases/download
# ENV GITHUB_APOC_URI=https://github.com/neo4j-contrib/neo4j-apoc-procedures/releases/download
# ENV NEO4J_GDS_URI=${GITHUB_GDS_URI}/${GDS_LIB_VERSION}/neo4j-graph-data-science-${GDS_LIB_VERSION}-standalone.jar
# ENV NEO4J_APOC_URI=${GITHUB_APOC_URI}/${APOC_LIB_VERSION}/apoc-${APOC_LIB_VERSION}-all.jar

# # Download Neo4j APOC libraries
# RUN sh -c 'cd /var/lib/neo4j/plugins && \
#     echo "Downloading APOC libraries..." \
#     curl -C- \
#         --location ${NEO4J_APOC_URI} \
#         --output $NEO4J_DIR/plugins/apoc-${APOC_LIB_VERSION}-all.jar'

# # Download Neo4j GDS libraries
# RUN sh -c 'cd /var/lib/neo4j/plugins && \
#     echo "Downloading Neo4j Graph Data Science libraries..." \
#     curl -C- \
#         --location ${NEO4J_GDS_URI} \
#         --output $NEO4J_DIR/plugins/neo4j-graph-data-science-${GDS_LIB_VERSION}-standalone.jar'

WORKDIR /var/lib/neo4j

# Make the backups directory
RUN mkdir -p /var/lib/neo4j/backups


# Mount the directories directly into Neo4j import directory
# VOLUME ../data/csv/ /var/lib/neo4j/import/
VOLUME /neo4j/data/ /var/lib/neo4j/data/
VOLUME /neo4j/conf/ /var/lib/neo4j/conf/
VOLUME /neo4j/plugins/ /var/lib/neo4j/plugins/
VOLUME /neo4j/logs/ /var/lib/neo4j/logs/
VOLUME /neo4j/backups/ /var/lib/neo4j/backups/

EXPOSE 7474 7473 7687
