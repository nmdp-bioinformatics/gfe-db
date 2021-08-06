FROM python:3.8

RUN apt update && \
    apt-get install -y \
        bc \
    && pip3 install --upgrade pip \
    && apt-get clean

RUN pip3 --no-cache-dir install --upgrade awscli

ENV GFE_BUCKET=gfe-db-4498

WORKDIR /gfe-db
RUN mkdir data

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

COPY . .

VOLUME data/ data/

CMD ["bash", "scripts/build.sh", "10"]