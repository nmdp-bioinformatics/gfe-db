FROM nmdpbioinformatics/gfe-base as gfe-graph-builder

WORKDIR /opt

ARG IMGT="3360"
ARG K=False
ARG AN=False

ENV RELEASES $IMGT
ENV KIR $K
ENV ALIGN $AN

RUN apk add curl

# Copy the build scripts to /opt
COPY bin/get_alignments.sh /opt/
COPY bin/build_gfedb.py /opt/build_gfedb.py
COPY bin/build.sh /opt/

RUN sh /opt/build.sh /opt


FROM neo4j:3.1 as neo4j-db-builder
COPY --from=gfe-graph-builder --chown=neo4j:neo4j /data/csv /csv

COPY bin/load_graph_docker.sh /opt/
RUN /opt/load_graph_docker.sh /csv

#
#

FROM neo4j:3.1
COPY --from=neo4j-db-builder --chown=neo4j:neo4j /var/lib/neo4j/gfedb/graph.db /data/databases/graph.db
