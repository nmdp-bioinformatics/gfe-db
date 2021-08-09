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
    - [Build Service](#build-service)
    - [Load Service](#load-service)
    - [Database Service](#database-service)
    - [CloudFormation Templates](#cloudformation-templates)
  - [Installation](#installation)
    - [Prerequisites](#prerequisites)
    - [Creating a Python Virtual Environment](#creating-a-python-virtual-environment)
    - [Environment Variables](#environment-variables)
- [Build and run Docker locally](#build-and-run-docker-locally)
    - [Memory Management](#memory-management)
  - [Deployment](#deployment)
  - [Clean Up](#clean-up)
  - [Authors](#authors)
  - [References & Links](#references--links)


## Project Structure
```
.
├── LICENSE
├── README.md
├── build                                 # Build service
│   ├── Dockerfile                        
│   ├── requirements.txt
│   ├── run.sh
│   ├── scripts
│   │   └── get_alignments.sh
│   └── src
│       ├── build_gfedb.py
│       └── constants.py
├── cfn                                   # CloudFormation templates
│   ├── database.yml
│   └── update-pipeline.yaml
├── load                                  # Load service
│   ├── Dockerfile
│   ├── cypher
│   │   ├── create_index.cyp
│   │   ├── delete_db.cyp
│   │   └── load.cyp
│   ├── requirements.txt
│   ├── run.sh
│   └── src
│       └── load_gfedb.py
├── neo4j                                 # Database
│   ├── Dockerfile
│   └── plugins
│       ├── apoc-4.2.0.2-all.jar
│       └── neo4j-graph-data-science-1.5.1-standalone.jar
└── notebooks                             # Development notebooks
    ├── 1.0-load-gfe-db.ipynb
    └── 1.0-refactor-gfedb_utils.ipynb
```

## Description
The `gfe-db` represents IPD-IMGT/HLA sequence data as GFE nodes and relationships in a Neo4j graph database. The architecture to run and update `gfe-db` contains 3 basic components:
- Build Service
- Load Service
- Database Service

### Build Service
The build service is triggered when a new IMGT/HLA version is released. AWS Batch is used to deploy a container to an EC2 instance which will run the build script and generate a dataset of CSVs. These are uploaded to S3 where they can be accessed by the load service. This service is located inside the `build/` directory.

### Load Service
The load service runs once the build service completes. For each CSV file generated by the build service, a pre-signed URL is created and inserted into the `LOAD CSV FROM ...` statement within the Cypher script. Each statement in the script is sent to the Neo4j server using the HTTP API. The load service runs until Neo4j is done loading. This service is located inside the `load/` directory.

### Database Service
Neo4j is deployed within a Docker container to an EC2 instance. Indexes and constraints are set to expedite transactions. The browser can be accessed at port 7474 of the public DNS server found in the EC2 console. This service is located inside the `neo4j/` directory.

### CloudFormation Templates
CloudFormation templates define the architecture that is deployed to AWS. The basic resources include:
- IAM permissions
- AWS Batch job definitions, queues and compute environments for build and load services
- StepFunctions state machine to orchestrate the build and load service
- ECR repositories to host the container images used for the build and load services
- EC2 Launch Template for deploying Neo4j

## Installation
Follow these steps to work with the files locally.
- Create separate virtual environments inside the `build/` and `load/` directories and (optional) create kernels to use each environment inside Jupyter Notebook.

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
To use the virtual environment inside Jupyter Notebook it is necessary to create a kernel.
```bash
# Install ipykernal
pip install ipykernel

# Add the kernel
python3 -m ipykernel install --user --name=<environment name>
```

### Environment Variables
Add a `.env` file with the following variables.
```bash
GFE_BUCKET=<value>
RELEASES=<value>
ALIGN=<value>
KIR=<value>
MEM_PROFILE=<value>
LIMIT=<value>
NEO4J_HOST=<value>
NEO4J_USERNAME=<value>
NEO4J_PASSWORD=<value>
```
Run the command to export environment variables to the shell.
```bash
set -a
source .env
set +a
```

*Important:* *Always use a `.env` file or AWS SSM Parameter Store for sensitive variables like credentials and API keys. Never hard-code them, including when developing. AWS will quarantine an account if any credentials get accidentally exposed and this will cause problems. **MAKE SURE `.env` IS LISTED IN `.gitignore`.**

<!-- ## Usage
Follow these steps in sequence to build and load `gfe-db` locally. Make sure that your environment variables are set correctly before proceeding. -->

<!-- ### Run Neo4j Docker
Build the Docker image as defined in the Dockerfile. See [Configuring Neo4j in Dockerfile](#Configuring-Neo4j-in-Dockerfile) for important configuration settings.
```
cd neo4j
docker build --tag gfe-db .
```
Run the container to start Neo4j in Docker.
```
# Run container to start Neo4j
docker run -d --name gfe-db \
  -v "$(pwd)"/../data/csv/:/var/lib/neo4j/import \
  -v "$(pwd)"/../neo4j/plugins:/var/lib/neo4j/plugins \
  -v "$(pwd)"/../neo4j/logs:/var/lib/neo4j/logs \
  -p 7474:7474 -p 7473:7473 \
  -p 7687:7687 gfe-db
```
Access the container logs during startup to check the status of Neo4j.
```bash
docker logs -f gfe-db
```
Stop and restart as needed.
```bash
# Stop container
docker stop gfe-db

# Start container
docker start gfe-db
``` -->

<!-- ### Build GFE dataset
Run the command to build the container for the build service.
```
# Build and run Docker locally
cd build
docker build --tag gfe-db-build-service .
```
Run the command to start the build.
```
docker run \
  --rm \
  -v "$(pwd)"/../data:/opt/data \
  -v "$(pwd)"/logs:/opt/app/logs \
  -e GFE_BUCKET='gfe-db-4498' \
  -e RELEASES='3440' \
  -e ALIGN='False' \
  -e KIR='False' \
  -e MEM_PROFILE='True' \
  -e LIMIT='100' \
  --name gfe-db-build-service \
  gfe-db-build-service:latest
``` -->
<!-- ```
# Load from Docker locally
cd load
docker build -t gfe-db-load-service .
docker run gfe-db-load-service:latest
``` -->

<!-- ### Load the dataset into Neo4j
Once the container is running, the Neo4j server is up, and the dataset has been created, run the command to load it into Neo4j.
```
bash bin/load_db.sh
``` -->

<!-- ### Running the GFE database in Neo4j 4.2 using Docker
This README outlines the steps for building and running a development version of `gfe-db` in a local Docker container. Docker will deploy an instance of Neo4j 4.2 including the [APOC](https://neo4j.com/labs/apoc/4.1/) and [Graph Data Science](https://neo4j.com/docs/graph-data-science/current/) plugins. GFE data is stored in the `data/csv/` directory which is mounted as an external volume within the container when run. This keeps the data outside the container so that it can be updated easily. -->

## Notebooks

### `1.0-load-gfe-db`
Python notebook for developing the load service using the Neo4j HTTP API, Requests and Boto3 libraries.

### `1.0-refactor-gfedb_utils`
Development notebook for refactoring `gfe-db` and the `build/src/build_gfedb.py` module used for building the CSV datasets.

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
`gfe-db` is deployed using Docker to an EC2 instance. Automated builds and loading of `gfe-db` on AWS is orchestrated using AWS Batch and StepFunctions. The infrastructure is defined using CloudFormation templates.

1. Make sure to update your AWS credentials in `~/.aws/credentials`
2. Create an S3 bucket and add the name to the `gfeBucket` parameter in `update-pipeline.yaml`.
3. Deploy the CloudFormation stacks.
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
```
2. Follow the instructions in each ECR repo to push the images to that repective repo.
3. Trigger an update using StepFunctions by starting an execution with the following input:
   ```json
   {
     "params": {
       "environment": {
         "RELEASES": "3450",
         "ALIGN": "False",
         "KIR": "False",
         "MEM_PROFILE": "False",
         "LIMIT": ""
       }
     }
   }
   ```
  Update the parameters to whatever is desired. Leaving `LIMIT` blank will build the entire GFE dataset (~30,000 alleles).


## Clean Up

To delete a stack and it's resources, run the command. S3 buckets and ECR repositories must be empty before they can be deleted.
```bash
aws cloudformation delete-stack --stack-name <stack name>
```

<!-- ### Local Clean-up
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
``` -->

## Troubleshooting
* Check your AWS credentials in `~/.aws/credentials`
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