FROM python:3.8-slim-buster

RUN apt update && \
    apt-get install -y \
        bc \
        curl \
        jq \
    && pip3 install --upgrade pip \
    && apt-get clean

WORKDIR /opt/app

COPY requirements.txt requirements.txt
RUN pip3 --no-cache-dir install -r requirements.txt

COPY . .

ENTRYPOINT ["bash", "run.sh"]