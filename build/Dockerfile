FROM python:3.8-slim-buster

RUN apt update && \
    apt-get install -y \
        bc \
        curl \
    && pip3 install --upgrade pip \
    && apt-get clean

RUN pip3 --no-cache-dir install --upgrade awscli

WORKDIR /opt/app
RUN mkdir /opt/data

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

COPY . .

VOLUME data/ /opt/data

ENTRYPOINT ["bash", "run.sh"]