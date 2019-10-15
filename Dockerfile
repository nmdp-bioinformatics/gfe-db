FROM gfe-base:latest as gfe-graph-builder

WORKDIR /opt

ENV NEO4J_HOME /var/lib/neo4j
ENV NEO4J_BIN /var/lib/neo4j/bin

ARG IMGT="3360"
ARG K=False
ARG AN=False

ENV RELEASES $IMGT
ENV KIR $K
ENV ALIGN $AN

COPY mod-imgt /mod-imgt/
COPY requirements.txt /opt/

RUN pip3 install -r /opt/requirements.txt

COPY bin/get_alignments.sh /opt/
COPY bin/build_gfedb.py /opt/build_gfedb.py
COPY bin/build.sh /opt/

COPY biopython /opt/biopython
RUN pip uninstall -y biopython && cd /opt/biopython/ && pip install .

RUN sh /opt/build.sh /opt


#
#EXPOSE 7474 7473 7687
#
#CMD sh /opt/load_graph.sh /opt


