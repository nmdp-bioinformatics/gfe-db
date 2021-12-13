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
  - [Services](#services)
    - [Infrastructure](#infrastructure)
    - [Database](#database)
    - [Pipeline](#pipeline)
  - [Installation](#installation)
    - [Prerequisites](#prerequisites)
    - [AWS Configuration](#aws-configuration)
    - [Environment Vairables](#environment-vairables)
  - [Deployment](#deployment)
  - [Local Development](#local-development)
    - [Creating a Python Virtual Environment](#creating-a-python-virtual-environment)
    - [Docker](#docker)
  - [Troubleshooting](#troubleshooting)
  - [Authors](#authors)
  - [References & Links](#references--links)

## Project Structure
```bash
.
├── LICENSE
├── Makefile
├── README.md
└── gfe-db
    ├── database             # Database service
    │   ├── Makefile
    │   ├── neo4j
    │   └── template.yaml
    ├── infrastructure       # Network infrastructure including VPC, SSM Parameters and Secrets
    │   ├── Makefile
    │   └── template.yaml
    └── pipeline             # Update pipeline including Batch jobs, StepFunctions, trigger
        ├── Makefile
        ├── config
        ├── functions
        │   ├── Makefile
        │   └── trigger
        ├── jobs
        │   ├── Makefile
        │   ├── build
        │   └── load
        └── template.yaml
```

## Description
The `gfe-db` represents IPD-IMGT/HLA sequence data as GFE nodes and relationships in a Neo4j graph database. Running this application will setup the following services in AWS:
- VPC and subnet
- Neo4j database server
- Update pipeline and trigger

## Services
The project organizes its resources by service so that deployments are decoupled (acheived using Makefiles). Shared configurations leverage SSM Parameter Store and Secrets Manager.

### Infrastructure
The infrastructure service deploys a VPC, public subnet, and common SSM parameters and secrets for the other services to use.

### Database
The database service deploys an EC2 instance hosting a Neo4j Docker container into a public subnet so that it can be accessed through a browser.

### Pipeline
The pipeline service automates updates of the database using a scheduled Lambda which will trigger a build and load of new data when it becomes available. The trigger Lambda watches the source data repository and triggers the pipeline when a new IMGT/HLA version is released. The pipeline uses a StepFunctions state machine to orchestrate the build and load steps using AWS Batch.

## Installation
Follow the steps to set the deployment environment.

### Prerequisites
* Python 3.8
* GNU Make 3.81
* AWS CLI
* SAM CLI
* Docker
* jq

### AWS Configuration
Valid AWS credentials must be available to AWS CLI and SAM CLI. The easiest way to do this is running `aws configure`, or by adding them to `~/.aws/credentials` and exporting the `AWS_PROFILE` variable to the environment.

For more information visit the documentation page:
[Configuration and credential file settings](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-files.html)

### Environment Vairables
Non-sensistive environment variables are handled by the root Makefile. Sensitive environment variables containing secrets like passwords and API keys must be exported to the environment first.

Create a `.env` file in the project root.
```bash
NEO4J_USERNAME=<value>
NEO4J_PASSWORD=<value>
GITHUB_PERSONAL_ACCESS_TOKEN=<value>
```

Source the variables to the environment.
```bash
set -a
source .env
set +a
```
*Important:* *Always use a `.env` file or AWS SSM Parameter Store or Secrets Manager for sensitive variables like credentials and API keys. Never hard-code them, including when developing. AWS will quarantine an account if any credentials get accidentally exposed and this will cause problems. **MAKE SURE `.env` IS LISTED IN `.gitignore`.**

## Deployment
Once an AWS profile is configured and environment variables are exported, the application can be deployed using `make`.
```bash
make deploy
```
It is also possible to deploy or update individual services.
```bash
# Deploy/update only the infrastructure service
make deploy.infrastructure

# Deploy/update only the database service
make deploy.database

# Deploy/update only the pipeline service
make deploy.pipeline
```
Note: It is recommended to only deploy from the project root. This is because common parameters are passed from the root Makefile to nested Makefiles. If a stack has not been changed, the deployment script will continue until it reaches a stack with changes and deploy that.

## Local Development

### Creating a Python Virtual Environment
When developing locally, you will need to create an individual virtual environment to run scripts in the `jobs` or `functions` directories, since they require different dependencies.
```bash
cd <specific job or function directory>
python3 -m venv .venv
source .venv/bin/activate
pip install -U pip
pip install -r requirements.txt
```

To use the virtual environment inside a Jupyter Notebook, first activate the virtual environment, then create a kernel for it.
```bash
# Install ipykernal
pip install ipykernel

# Add the kernel
python3 -m ipykernel install --user --name=<environment name>

# Remove the kernel
jupyter kernelspec uninstall <environment name>
```

### Docker
Build the Docker image as defined in the Dockerfile.
```bas
cd <directory>
docker build --tag gfe-db .
```
Run the container to start Docker. Specify volumes, ports, or environment variables if necessary.
```
# Run Neo4j Docker
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

# Access the container's shell
docker exec -it <container> bash
```

<!-- ### Build GFE dataset
Run the command to build the container for the build service.
```
# Build and run Docker locally
cd build
docker build --tag gfe-db-build-service .
```
Run the command to start the build. (Requires an S3 bucket)
```
docker run \
  --rm \
  -v "$(pwd)"/../data:/opt/data \
  -v "$(pwd)"/logs:/opt/app/logs \
  -e GFE_BUCKET='<S3 bucket name>' \
  -e RELEASES='3440' \
  -e ALIGN='False' \
  -e KIR='False' \
  -e MEM_PROFILE='True' \
  -e LIMIT='100' \
  --name gfe-db-build-service \
  gfe-db-build-service:latest
``` 
```
# Load from Docker locally
cd load
docker build -t gfe-db-load-service .
docker run gfe-db-load-service:latest
```

### Load the dataset into Neo4j
Once the container is running, the Neo4j server is up, and the dataset has been created, run the command to load it into Neo4j.
```
bash bin/load_db.sh
``` -->

<!-- ### Running the GFE database in Neo4j 4.2 using Docker
This README outlines the steps for building and running a development version of `gfe-db` in a local Docker container. Docker will deploy an instance of Neo4j 4.2 including the [APOC](https://neo4j.com/labs/apoc/4.1/) and [Graph Data Science](https://neo4j.com/docs/graph-data-science/current/) plugins. GFE data is stored in the `data/csv/` directory which is mounted as an external volume within the container when run. This keeps the data outside the container so that it can be updated easily. -->



<!-- ### `1.0-load-gfe-db`
Python notebook for developing the load service using the Neo4j HTTP API, Requests and Boto3 libraries.

### `1.0-refactor-gfedb_utils`
Development notebook for refactoring `gfe-db` and the `build/src/build_gfedb.py` module used for building the CSV datasets. -->

<!-- ## Running Tests -->

<!-- ## Configuring Neo4j in Dockerfile
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
``` -->

<!-- ## Deployment
`gfe-db` is deployed using Docker to an EC2 instance. Automated builds and loading of `gfe-db` on AWS is orchestrated using AWS Batch and StepFunctions. The infrastructure is defined using CloudFormation templates.

1. Make sure to update your AWS credentials in `~/.aws/credentials`.
2. Deploy the CloudFormation stacks and set the default region to `us-east-1` (there are issues with S3 pre-signed URLs in the other regions).
   ```bash
   cd gfe-db
   bash deploy.sh
   ```
3. In the AWS ECR console, follow the instructions in each ECR repo to build, tag and push the images to that repo.
4. In the Neo4j browser, run the `load/cypher/create_index.cyp` script.
5. Trigger an update using StepFunctions by starting an execution with the following input:
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
6. Get the `Neo4jDatabaseEndpoint` parameter from the `database-stack.yml` and load it in the browser on port 7474:
   ```bash
   # Example
   18.215.230.187:7474
   ```
   The graph should be loaded once the StepFunctions completes.


## Clean Up
To delete a stack and it's resources, use the CloudFormation console or run the command. S3 buckets and ECR repositories must be empty before they can be deleted.
```bash
aws cloudformation wait stack-delete-complete --stack-name <stack name>
``` -->

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