FROM python:3.8

RUN apt update && \
    apt-get install -y \
        bc \
    && pip3 install --upgrade pip \
    && apt-get clean

# Overidden by AWS Batch
ENV GFE_BUCKET=gfe-db-4498
ENV RELEASES="3440"
ENV ALIGN=True
ENV KIR=False
ENV MEM_PROFILE=True

WORKDIR /gfe-db

COPY requirements.txt requirements.txt
RUN pip3 --no-cache-dir install -r requirements.txt

COPY . .

ENTRYPOINT ["bash", "load_db.sh"]