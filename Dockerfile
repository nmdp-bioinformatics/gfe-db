FROM neo4j:4.2

# RUN apt-get update \
#     && apt-get install -y curl openssl apt-utils zip unzip

ENV NEO4J_AUTH=neo4j/gfedb \
    NEO4J_ACCEPT_LICENSE_AGREEMENT=yes \
    NEO4J_dbms_memory_heap_initial__size=2G \
    NEO4J_dbms_memory_heap_max__size=2G \
    NEO4J_dbms_memory_pagecache_size=1G \
    NEO4J_dbms_security_procedures_unrestricted=apoc.*,gds.* \
    NEO4J_dbms_security_allow__csv__import__from__file__urls=true
    # NEO4J_apoc_import_file_enabled=true \
    # NEO4J_apoc_import_file_use__neo4j__config=true \
    # NEO4J_apoc_export_file_enabled=true 

# ENV APOC_URI https://github.com/neo4j-contrib/neo4j-apoc-procedures/releases/download/${APOC_VERSION}/apoc-${APOC_VERSION}-all.jar
# RUN sh -c 'cd /var/lib/neo4j/plugins && curl -L -O "${APOC_URI}"'

# TO DO: add the downloaded plugins (APOC, GDS)
COPY data/csv/ /var/lib/neo4j/import/

EXPOSE 7474 7473 7687

CMD ["neo4j"]
