FROM 682793961433.dkr.ecr.us-east-1.amazonaws.com/gfedb

COPY bin/* opt/
COPY requirements.txt opt/

WORKDIR /opt

ENV NEO4J_HOME /var/lib/neo4j
ENV NEO4J_BIN /var/lib/neo4j/bin
ENV NEO4J_CONF /opt/conf

ARG IMGT="3170,3190,3200,3210,3240,3250,3310"
ARG K=False
ARG AN=False

ENV RELEASES $IMGT
ENV KIR $K
ENV ALIGN $AN

RUN pip3 install -r /opt/requirements.txt

RUN sh /opt/build.sh /opt

EXPOSE 7474 7473 7687

CMD sh /opt/load_graph.sh /opt


