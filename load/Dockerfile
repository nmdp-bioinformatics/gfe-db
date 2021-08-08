FROM python:3.8

RUN apt update && \
    apt-get install -y \
        bc \
    && pip3 install --upgrade pip \
    && apt-get clean

# Overidden by AWS Batch
ENV GFE_BUCKET=gfe-db-4498
ENV RELEASES=3440
ENV NEO4J_HOST=44.192.54.30
ENV NEO4J_USERNAME=neo4j
ENV NEO4J_PASSWORD=gfedb

WORKDIR /opt/app

COPY requirements.txt requirements.txt
RUN pip3 --no-cache-dir install -r requirements.txt

COPY . .

ENTRYPOINT ["bash", "run.sh"]