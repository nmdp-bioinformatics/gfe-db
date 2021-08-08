FROM python:3.8

RUN apt update && \
    apt-get install -y \
        bc \
    && pip3 install --upgrade pip \
    && apt-get clean

RUN pip3 --no-cache-dir install --upgrade awscli

ENV AWS_ACCESS_KEY_ID=''
ENV AWS_SECRET_ACCESS_KEY=''

# Overidden by AWS Batch
ENV RELEASES=3440
ENV ALIGN=True
ENV KIR=False
ENV MEM_PROFILE=True

WORKDIR /gfe-db
RUN mkdir data

COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt

COPY . .

VOLUME data/ data/

ENTRYPOINT ["bash", "scripts/build.sh"]