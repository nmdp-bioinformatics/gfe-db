FROM mhalagan1nmdp/gfe-base:latest 

COPY bin/*.* /opt/
COPY mod-imgt /mod-imgt/
COPY requirements.txt /opt/

WORKDIR /opt

ENV NEO4J_HOME /var/lib/neo4j
ENV NEO4J_BIN /var/lib/neo4j/bin

ARG IMGT="3360"
ARG K=False
ARG AN=False

ENV RELEASES $IMGT
ENV KIR $K
ENV ALIGN $AN

RUN pip3 install -r /opt/requirements.txt

RUN sh /opt/build.sh /opt

EXPOSE 7474 7473 7687

CMD sh /opt/load_graph.sh /opt


