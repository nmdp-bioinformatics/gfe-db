FROM neo4j:5.15-community-bullseye

ENV NEO4J_AUTH=none
ENV NEO4J_PLUGINS='["apoc", "apoc-extended", "graph-data-science"]'

ARG NEO4J_DATABASE_NAME

RUN sed -i "s/#initial.dbms.default_database=neo4j/initial.dbms.default_database=$NEO4J_DATABASE_NAME/g" "${NEO4J_HOME}/conf/neo4j.conf"
RUN sed -i "s/#dbms.security.procedures.allowlist=apoc.coll.*,apoc.load.*,gds.*/dbms.security.procedures.allowlist=apoc.*,gds.*/g" "${NEO4J_HOME}/conf/neo4j.conf"


COPY neo4j/data "${NEO4J_HOME}/data"