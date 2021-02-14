FROM neo4j:4.2

# RUN apt-get update \
#     && apt-get install -y curl openssl apt-utils zip unzip

ENV NEO4J_AUTH=neo4j/gfedb \
    NEO4J_ACCEPT_LICENSE_AGREEMENT=yes
    # APOC_VERSION=4.1.0.6

# ENV APOC_URI https://github.com/neo4j-contrib/neo4j-apoc-procedures/releases/download/${APOC_VERSION}/apoc-${APOC_VERSION}-all.jar
# RUN sh -c 'cd /var/lib/neo4j/plugins && curl -L -O "${APOC_URI}"'

COPY data/csv/ /var/lib/neo4j/import/

EXPOSE 7474 7473 7687

CMD ["neo4j"]
