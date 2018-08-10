FROM openjdk:8-jre-alpine

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

RUN echo $RELEASES $KIR $ALIGN

RUN apk add --no-cache --quiet \
    bash \
    curl \
    perl

RUN apk update && apk add gcc g++ make patch perl-dev wget
RUN curl -L http://xrl.us/cpanm > /bin/cpanm && chmod +x /bin/cpanm
RUN cpanm App::cpm

RUN cpanm install --force \
  Bio::Perl 

RUN apk add --no-cache python3 && \
    python3 -m ensurepip && \
    rm -r /usr/lib/python*/ensurepip && \
    pip3 install --upgrade pip setuptools && \
    if [ ! -e /usr/bin/pip ]; then ln -s pip3 /usr/bin/pip ; fi && \
    if [[ ! -e /usr/bin/python ]]; then ln -sf /usr/bin/python3 /usr/bin/python; fi && \
    rm -r /root/.cache

RUN apk add libxml2 libxslt libxml2-dev libxslt-dev

RUN apk add --no-cache python3-dev && \
    apk add --no-cache --virtual .build-deps g++ && \
    ln -s /usr/include/locale.h /usr/include/xlocale.h

RUN pip3 install -r /opt/requirements.txt

ENV PERL5LIB=/usr/local/lib/perl5
ENV PATH=/usr/local/bin:$PATH

ENV NEO4J_SHA256 f0d79b4a98672dc527b708113644b8961ba824668c354e61dc4d2a16d8484880
ENV NEO4J_TARBALL neo4j-community-3.1.3-unix.tar.gz
ARG NEO4J_URI=http://dist.neo4j.org/neo4j-community-3.1.3-unix.tar.gz

RUN curl --fail --silent --show-error --location --remote-name ${NEO4J_URI} \
    && echo "${NEO4J_SHA256}  ${NEO4J_TARBALL}" | sha256sum -csw - \
    && tar --extract --file ${NEO4J_TARBALL} --directory /var/lib \
    && mv /var/lib/neo4j-* /var/lib/neo4j \
    && rm ${NEO4J_TARBALL}

RUN sh /opt/build.sh /opt

EXPOSE 7474 7473 7687

CMD sh /opt/load_graph.sh /opt


