FROM python:3.8

RUN apt update && \
    apt-get install -y \
        bc \
    && pip3 install --upgrade pip \
    && apt-get clean

RUN pip3 --no-cache-dir install --upgrade awscli

WORKDIR /gfe-db
RUN mkdir data

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

COPY . .

VOLUME data/ data/

ENTRYPOINT ["bash", "run.sh"]