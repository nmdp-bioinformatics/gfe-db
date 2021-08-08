gfe-db
======

Graph database representing IPD-IMGT/HLA sequence data as GFE.

<!-- Need to update image -->
<!-- ![GFE Schema](img/gfe-schema.png) -->

## Table of Contents
- [gfe-db](#gfe-db)
  - [Table of Contents](#table-of-contents)
  - [Project Structure](#project-structure)
  - [Description](#description)
    - [New Features](#new-features)
  - [Installation](#installation)
    - [Prerequisites](#prerequisites)
    - [Creating a Python Virtual Environment](#creating-a-python-virtual-environment)
    - [Environment Variables](#environment-variables)
  - [Usage](#usage)
    - [Run Neo4j Docker](#run-neo4j-docker)
    - [Build GFE dataset](#build-gfe-dataset)
  - [Configuring Neo4j in Dockerfile](#configuring-neo4j-in-dockerfile)
    - [Username & Password](#username--password)
    - [Memory Management](#memory-management)
  - [Deployment](#deployment)
- [Deploy update pipeline](#deploy-update-pipeline)
  - [Clean Up](#clean-up)
    - [Local Clean-up](#local-clean-up)
  - [Troubleshooting](#troubleshooting)
  - [Authors](#authors)
  - [References & Links](#references--links)


## Project Structure
```
├── bin
│   ├── build.sh                    # Executable scripts
│   ├── get_alignments.sh           # Alignments are included by default
│   ├── load_db.sh                  # Loads multiple IMGT/HLA versions
│   └── set_env.sh                  # Exports environment variables
├── (data)                          # Created during build step
│   ├── 3360                        # Alignments
│   ├── csv                         # CSVs loaded into Neo4j for each IMGT release
│   └── hla.3360.dat                # Allele data                
├── neo4j                           # Neo4j load scripts
│   ├── create_index.cyp
│   ├── delete_db.cyp
│   ├── load.cyp                    # Merges new nodes with existing
│   ├── logs
│   └── plugins
├── notebooks                       # Development jupyter notebooks
├── src                             # Build modules
│   ├── __init__.py
│   ├── constants.py
│   ├── gfedb.py                    
│   └── gfedb_utils.py              
├── .dockerignore                   # Files for Docker to ignore
├── .gitignore                      # Files for git to ignore
├── Dockerfile                      # Docker image for Neo4j 4.2
├── LICENSE
├── README.md                       # Instructions for this workflow
├── EC2INSTRUCTIONS.md              # Instructions for deployment on EC2
└── requirements.txt                # Python dependencies
```

## Description
The `gfe-db` represents IPD-IMGT/HLA sequence data as GFE nodes and relationships in a Neo4j graph database.

### New Features
* This version uses a Neo4j 4.2 Docker image.
* Multiple IMGT/HLA releases can be loaded into the same graph.
* Schema is updated.

Please feel free to open issues regarding specific bugs and feature requests.

## Installation

### Prerequisites
* Python 3.8
* Docker
* AWS CLI

### Creating a Python Virtual Environment
Run these commands to create a virtual environment that will install the libraries listed in `requirements.txt`.
```bash
# Create .venv and activate
python3 -m venv .venv
source .venv/bin/activate

# Update pip to avoid conflicts
pip install -U pip

# Install libraries
pip install -r requirements.txt
```

### Environment Variables
The build script requires specific variables to be present in the environment. To set these, follow these steps.
1. Create a `scripts/set_env.sh` file in the project root
2. Add these variables to the `set_env.sh` file
```bash
#!/bin/bash

export GFE_BUCKET=<S3 bucket name>
export RELEASES="<release version>"
export ALIGN=<boolean>
export KIR=<boolean>
export MEM_PROFILE=<boolean>
```
3. Run the command.
```bash
source scripts/set_env.sh
```

*Important:* *Always use a `.env` file or AWS SSM Parameter Store for sensitive variables like credentials and API keys. Never hard-code them, including when developing. AWS will quarantine an account if any credentials get accidentally exposed and this will cause problems. **MAKE SURE `.env` IS LISTED IN `.gitignore`.**

## Usage
Follow these steps in sequence to build and load `gfe-db` locally.

### Run Neo4j Docker
Build the Docker image as defined in the Dockerfile. See [Configuring Neo4j in Dockerfile](#Configuring-Neo4j-in-Dockerfile) for important configuration settings.
```
docker build --tag gfe-db .
```
Run the container to start Neo4j in Docker.
```
# Run container to start Neo4j
docker run -d --name gfe \
  -v "$(pwd)"/data/csv/:/var/lib/neo4j/import \
  -v "$(pwd)"/neo4j/plugins:/var/lib/neo4j/plugins \
  -v "$(pwd)"/neo4j/logs:/var/lib/neo4j/logs \
  -p 7474:7474 -p 7473:7473 \
  -p 7687:7687 gfe-db
```
If desired, access the container logs during startup. This will indicate when Neo4j is ready.
```bash
docker logs -f gfe
```
Stop and restart as needed.
```bash
# Stop container
docker stop gfe

# Start container
docker start gfe
```

### Build GFE dataset
To build gfe-db locally, run the command.
```bash
bash bin/build.sh

# Build 100 alleles (the full database contains around 30,000 alleles)
bash bin/build.sh 100
```

<!-- ```
# Build and run Docker locally
docker build -t gfe-db-build-service build/
docker run -v "$(pwd)"/data:/data gfe-db-build-service:latest
``` -->

### Load the dataset into Neo4j
Once the container is running, the Neo4j server is up, and the dataset has been created, run the command to load it into Neo4j.
```
bash bin/load_db.sh
```

<!-- ### Running the GFE database in Neo4j 4.2 using Docker
This README outlines the steps for building and running a development version of `gfe-db` in a local Docker container. Docker will deploy an instance of Neo4j 4.2 including the [APOC](https://neo4j.com/labs/apoc/4.1/) and [Graph Data Science](https://neo4j.com/docs/graph-data-science/current/) plugins. GFE data is stored in the `data/csv/` directory which is mounted as an external volume within the container when run. This keeps the data outside the container so that it can be updated easily. -->

## Notebooks

### `0.0-load-gfe-db`
Contains code to load the Neo4j database using Python and the Requests module.

### `1.0-refactor-gfedb_utils`
Development notebook for refactoring `gfe-db` and the `gfe-db_utils.py` module used for building the data.

### Adding a kernel spec to Jupyter Notebook
To use the virtual environment inside Jupyter Notebook it is necessary to create a kernel.
```bash
# Install ipykernal
pip install ipykernel

# Add the kernel
python3 -m ipykernel install \
  --user \
  --name gfe-db \
  --display-name "gfe-db"

# If running the notebooks, install these libraries not used by the application (optional)
pip install python-dotenv
```

To remove the kernel spec, run the command.
```bash
jupyter kernelspec uninstall gfe-db
```

<!-- ## Running Tests -->

## Configuring Neo4j in Dockerfile
Configuration settings for Neo4j are passed through environment variables in the Dockerfile.
### Username & Password
The username and password is set as follows:
```Dockerfile
# Dockerfile
ENV NEO4J_AUTH=neo4j/gfedb
```
### Memory Management
Optimal memory for Neo4j depends on available RAM. Loading and querying a larger dataset will require more memory allocated. Make sure that the Docker daemon is configured to handle whatever values are given here. For more information on memory management in Neo4j, see the [Neo4j Operations Manual](https://neo4j.com/docs/operations-manual/current/performance/memory-configuration/).
```Dockerfile
# Dockerfile; Rebuild the image after updating these
ENV NEO4J_dbms_memory_heap_initial__size=2G
ENV NEO4J_dbms_memory_heap_max__size=2G
```

## Deployment
`gfe-db` can be deployed using Docker to an EC2 instance. Automated builds and loading of `gfe-db` on AWS is orchestrated using AWS Batch and StepFunctions. The infrastructure is defined using CloudFormation templates.

```bash
# Deploy database server
aws cloudformation deploy \
   --template-file cfn/database.yml \
   --stack-name gfe-db \
   --capabilities CAPABILITY_NAMED_IAM

# Deploy build service
aws cloudformation deploy \
  --template-file cfn/update-pipeline.yaml \
  --stack-name gfe-db-update-pipeline \
  --capabilities CAPABILITY_NAMED_IAM

# Deploy CI/CD
aws cloudformation deploy \
  --template-file cfn/cicd.yaml \
  --stack-name gfe-db-cicd
```
<!-- ```
# Deploy update pipeline
aws cloudformation deploy \
  --template-file cfn/pipeline.yaml \
  --stack-name gfe-db-pipeline \
  --capabilities CAPABILITY_NAMED_IAM
``` -->

To delete a stack and it's resources, run the command. S3 buckets must be empty before they can be deleted.
```bash
aws cloudformation delete-stack --stack-name <stack name>
```

## Clean Up

### Local Clean-up
Delete the Docker container.
```bash
docker stop gfe
docker rm gfe
```
Delete the Docker image.
```bash
# List the images and get the IMAGE IDs for gfe-db and neo4j
docker image ls

# Remove the images for gfe-db:latest and neo4j:4.2 by the IMAGE ID
docker image rm <IMAGE ID> <IMAGE ID>
```
To delete the Jupyter Notebook kernel, run the command.
```bash
jupyter kernelspec uninstall gfe-db
```

## Troubleshooting
* Check that the environment variables have been exported
* Sometimes Neo4j sets permissions on mounted volumes. To get around this run this from the project root:
  ```bash
  sudo chmod -R 777 .
  ```
* Check that the virtual environment is activated: `source .venv/bin/activate`
* Check that requirements are installed: `pip install -r requirements.txt`
* Check that Python 3.8 is being used

## Authors
**Primary Contact:** Gregory Lindsey ([@abk7777](https://github.com/abk7777))

## References & Links
 * [hub.docker.com/r/nmdpbioinformatics/service-gfe-submission](https://hub.docker.com/r/nmdpbioinformatics/service-gfe-submission)
 * [service-gfe-submission.readthedocs.io](https://service-gfe-submission.readthedocs.io/en/latest/index.html)
 * [github.com/nmdp-bioinformatics/service-feature](https://github.com/nmdp-bioinformatics/service-feature)
 * [github.com/nmdp-bioinformatics/HSA](https://github.com/nmdp-bioinformatics/HSA)
 * [bioinformatics.bethematchclinical.org](https://bioinformatics.bethematchclinical.org)
 * [feature.nmdp-bioinformatics.org](https://feature.nmdp-bioinformatics.org)
 * [gfe.b12x.org](http://gfe.b12x.org)

-----------------
<br>
<p align="center">
  <img src="https://bethematch.org/content/site/images/btm_logo.png">
</p>