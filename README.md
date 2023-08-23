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
  - [Architecture](#architecture)
    - [Base Infrastructure](#base-infrastructure)
    - [Database](#database)
    - [Data Pipeline](#data-pipeline)
  - [Deployment](#deployment)
    - [Quick Start](#quick-start)
    - [Prerequisites](#prerequisites)
    - [Environment](#environment)
      - [AWS Credentials](#aws-credentials)
      - [Shell Variables](#shell-variables)
    - [Makefile Usage](#makefile-usage)
      - [Makefile Command Reference](#makefile-command-reference)
  - [Managing Configuration](#managing-configuration)
    - [Database Configuration](#database-configuration)
      - [Neo4j](#neo4j)
      - [Shell Scripts](#shell-scripts)
      - [Cypher Scripts](#cypher-scripts)
      - [SSL Policy](#ssl-policy)
    - [Data Pipeline Config](#data-pipeline-config)
      - [Invocation Input](#invocation-input)
      - [IMGT/HLA Release Versions State](#imgthla-release-versions-state)
  - [Loading Neo4j](#loading-neo4j)
    - [Clean Up](#clean-up)
  - [Backup \& Restore](#backup--restore)
    - [Backups](#backups)
    - [Restore](#restore)
  - [Local Development](#local-development)
    - [Creating a Python Virtual Environment](#creating-a-python-virtual-environment)
  - [Documentation](#documentation)
    - [Editing and Building the Documentation](#editing-and-building-the-documentation)
  - [Troubleshooting](#troubleshooting)
  - [Authors](#authors)
  - [References \& Links](#references--links)

## Project Structure
```bash
.
├── LICENSE
├── CONTRIBUTING.md
├── Makefile                                        # Use the root Makefile to deploy, delete and manage resources and configuration
├── README.md
├── docs                                            # Sphinx documentation
├── (<stage>-gfe-db-<region>-neo4j-key.pem)         # EC2 key pair for SSH access to Neo4j server, created on deployment
├── requirements-dev.txt                                # Python requirements for local development
├── requirements-docs.txt                                    # Python requirements for documentation
└── gfe-db
    # Database layer
    ├── database                                    # Neo4j database server, configuration and automation
    │   ├── Makefile
    │   ├── amazon-cloudwatch-agent
    │   │   └── amazon-cloudwatch-agent.json        # Sends EC2 logs to CloudWatch Logs for monitoring
    │   ├── neo4j
    │   │   ├── cypher                              # Cypher scripts for initialization and loading
    │   │   └── neo4j.template                      # Neo4j server configuration file
    │   ├── scripts                                 # Shell scripts for automation, loading, backup & restore
    │   └── template.yaml
    # Base Infrastructure layer
    ├── infrastructure                               # VPC and subnets, S3 bucket, SSM Parameters and Secrets
    │   ├── Makefile
    │   └── template.yaml
    # Data Pipeline layer
    └── pipeline                                     # Data pipeline including Batch job, Lambda functions & state machine
        ├── Makefile
        ├── config                                   # JSON files for storing app state
        │   ├── IMGTHLA-repository-state.json        # IMGT/HLA release version state
        │   └── pipeline-input.json                  # Default pipeline input parameters for scheduled invocations
        ├── functions
        │   ├── environment.json                     # Lambda configurations
        │   ├── invoke_load_script                   # Invokes Run Command to download and execute a schell script on EC2
        │   └── invoke_pipeline                      # Invokes the state machine
        ├── jobs                                     # Build job, triggered when a new IMGT/HLA version is released
        │   ├── Makefile
        │   └── build
        ├── statemachines                            # State machine for loading data into Neo4j
        └── template.yaml
```

## Description
The `gfe-db` represents IPD-IMGT/HLA sequence data as GFE nodes and relationships in a Neo4j graph database. This application deploys and configures AWS resources for the GFE database and an automated data pipeline for updates.

<br>
<p align="center">
  <img src="docs/source/_static/img/schema-alpha-v230305.png" alt="gfe-db schema" height="75%" width="75%">
</p>

## Architecture
<br>
<p align="center">
  <img src="docs/source/_static/img/arch-v230305.png" alt="gfe-db architecture diagram">
</p>

`gfe-db` architecture is organized by 3 layers each with its own Makefile:
1) Base Infrastructure 
2) Database 
3) Data pipeline
 
This allows the database and pipeline layers to be decoupled from each other and deployed or destroyed independently without affecting the other. Common configuration parameters are shared between resources using environment variables, JSON files, AWS SSM Paramter Store and Secrets Manager.

### Base Infrastructure
The base infrastructure layer deploys a VPC (optional), public subnet (optional), S3 bucket, Elastic IP and common SSM parameters and secrets for the other services to use.

### Database
The database layer deploys an EC2 instance running the Bitnami Neo4j AMI (Ubuntu 18.04) into a public subnet. An A record is required for a pre-existing Route53 domain and hosted zone so that SSL can be used to connect to Neo4j. During database deployment the SSL certificate is created and Cypher queries are run to create constraints and indexes, which help speed up loading and ensure data integrity. Neo4j is ready to be accessed through a browser once the instance has booted sucessfully.

During loading, a Lambda function calls the SSM Run Command API to execute bash scripts on the database instance. These scripts communicate with the Step Functions API to retrieve the job parameters, fetch the CSVs from S3 and populate the graph in Neo4j.

It is also possible to backup & restore to and from S3 by specific date checkpoints.

### Data Pipeline
The data pipeline layer automates integration of newly released IMGT/HLA data into Neo4j using a scheduled Lambda which watches the source data repository and invokes the build and load processes when it detects a new IMGT/HLA version in the upstream repository. The pipeline consists of a Step Functions state machine which orchestrates the build and load stages. The build process employs a Batch jobs to generate an intermediate set of CSV files. The load process leverages SSM Run Command to copy the CSV files to the Neo4j server and execute Cypher statements directly on the server (server-side loading). When loading the full dataset of 35,000+ alleles, the build step will generally take around 15 minutes, however the load step can take an hour or more.

## Deployment
It is possible to deploy gfe-db within it's own VPC, or to connect it to an external VPC by specigying `CREATE_VPC=true/false`.

### Quick Start
These list outline the basic steps for deployments. For more details please see the following sections.

**Using external VPC**
1. Retrieve the VPC ID and subnet ID from the AWS console or using the AWS CLI.
2. Purchase or designate a domain in Route53 and create a hosted zone with an A record for the subdomain. You can use `0.0.0.0` for the A record because it will be updated later by the deployment script.
3. Acquire a subscription for the Bitnami Neo4j AMI through [AWS Marketplace](https://aws.amazon.com/marketplace/pp/prodview-v47qqrn2yy7ie?sr=0-4&ref_=beagle&applicationId=AWSMPContessa).
4. [Install prerequisites](#Prerequisites).
5. [Set environment variables](#environment) including the ones from the previous steps. You must store these in a file named `.env.<stage>`, for example `.env.dev` or `.env.prod`:
    - CREATE_VPC=false
    - VPC_ID
    - PUBLIC_SUBNET_ID
    - HOSTED_ZONE_ID
    - HOST_DOMAIN
    - SUBDOMAIN
    - NEO4J_AMI_ID
6. Check the [config JSONs](#data-pipeline-config) (parameters and state) and edit the values as desired.
7. Run `STAGE=<stage> make deploy` to deploy the stacks to AWS.
8. Run `STAGE=<stage> make database.load.run releases=<version>` to load the Neo4j, or `STAGE=<stage> make database.load.run releases=<version> limit=<limit>` to run with a limited number of alleles.
9. Run `STAGE=<stage> make database.get.credentials` to get the username and password for Neo4j.
10. Run `STAGE=<stage> make database.get.endpoint` to get the URL for Neo4j and navigate to the Neo4j browser at the subdomain and host domain, for example `https://gfe-db.cloudftl.com:7473/browser/`.

**Creating a new VPC**
1. Purchase or designate a domain in Route53 and create a hosted zone with an A record for the subdomain. You can use `0.0.0.0` for the A record because it will be updated later by the deployment script.
2. Acquire a subscription for the Bitnami Neo4j AMI through [AWS Marketplace](https://aws.amazon.com/marketplace/pp/prodview-v47qqrn2yy7ie?sr=0-4&ref_=beagle&applicationId=AWSMPContessa).
3. [Install prerequisites](#Prerequisites).
4. [Set environment variables](#environment) including the ones from the previous steps. You must store these in a file named `.env.<stage>`, for example `.env.dev` or `.env.prod`:
    - CREATE_VPC=true
    - HOSTED_ZONE_ID
    - HOST_DOMAIN
    - SUBDOMAIN
    - NEO4J_AMI_ID
5. Check the [config JSONs](#data-pipeline-config) (parameters and state) and edit the values as desired.
6. Run `STAGE=<stage> make deploy` to deploy the stacks to AWS.
7. Run `STAGE=<stage> make database.load.run releases=<version>` to load the Neo4j, or `STAGE=<stage> make database.load.run releases=<version> limit=<limit>` to run with a limited number of alleles.
8. Run `STAGE=<stage> make database.get-credentials` to get the username and password for Neo4j.
9. Run `STAGE=<stage> make database.get.endpoint` to get the URL for Neo4j and navigate to the Neo4j browser at the subdomain and host domain, for example `https://gfe-db.cloudftl.com:7473/browser/`.

### Prerequisites
Please refer to the respective documentation for specific installation instructions.
* Route53 domain, hosted zone, and A record
* VPC & Public Subnet (if using external VPC)
* Bitnami Neo4j AMI subscription and AMI ID
* GNU Make 3.81
* coreutils (optional but recommended)
* AWS CLI
* SAM CLI
* Docker
* jq
* Python 3.9+ (if developing locally)

### Environment

####  AWS Credentials
Valid AWS credentials must be available to AWS CLI and SAM CLI. The easiest way to do this is with the following steps.
1. Run `aws configure` and follow the prompts, or copy/paste them into `~/.aws/credentials` 
2. Export the `AWS_PROFILE` variable for the chosen profile to the shell environment.
```bash
export AWS_PROFILE=default
```

For more information visit the documentation page:
[Configuration and credential file settings](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-files.html)

#### Shell Variables
These variables must be defined before running Make. The best way to set these variables is with a `.env.<stage>` file following this structure.
```bash
# .env.<stage>
AWS_PROFILE=<aws_profile> # Include profile if stacks are in a different accounts
STAGE=<dev or prod>
APP_NAME=<app_name>
AWS_REGION=<aws_region>
ADMIN_EMAIL=<email>
SUBSCRIBE_EMAILS=<email>,<email>,<email>,...
GITHUB_REPOSITORY_OWNER=<github_owner>
GITHUB_REPOSITORY_NAME=<github_repo_name>
HOST_DOMAIN=<fully_qualified_domain_name>
CREATE_VPC=<true or false>
VPC_ID=<vpc_id> # if CREATE_VPC=false
PUBLIC_SUBNET_ID=<public_subnet_id> # if CREATE_VPC=false
HOSTED_ZONE_ID=<hosted_zone_id>
SUBDOMAIN=<subdomain>
NEO4J_AMI_ID=<ami_id> # Requires AWS Marketplace subscription
APOC_VERSION=<apoc_version>
GDS_VERSION=<gds_version>
GITHUB_PERSONAL_ACCESS_TOKEN=<secret>
```

| Variable Name                | Example Value                      | Type   | Description                                      |
| ---------------------------- | ---------------------------------- | ------ | ------------------------------------------------ |
| AWS_PROFILE                  | <aws_profile>                      | string | AWS profile for deployment.                      |
| STAGE                        | dev                                | string | The stage of the application.                    |
| APP_NAME                     | gfe-db                             | string | The name of the application.                     |
| AWS_REGION                   | us-east-1                          | string | The AWS region to deploy to.                     |
| ADMIN_EMAIL                  | user@company.com                   | string | Admin's email required for SSL certificate.      |
| SUBSCRIBE_EMAILS             | user@company.com,user2@company.com | string | Comma-separated list of emails for notifications |
| GITHUB_REPOSITORY_OWNER      | <github_owner>                     | string | GitHub repository owner.                         |
| GITHUB_REPOSITORY_NAME       | <github_repo_name>                 | string | GitHub repository name.                          |
| HOST_DOMAIN                  | example.com                        | string | The domain to deploy to.                         |
| CREATE_VPC                   | true or false                      | string | Whether to create a new VPC.                     |
| HOSTED_ZONE_ID               | Z1234567890ABCDEF                  | string | The ID of the hosted zone to deploy to.          |
| SUBDOMAIN                    | gfe-db                             | string | The subdomain to deploy to.                      |
| NEO4J_AMI_ID                 | ami-0b9a2b6b1c5b8b5b9              | string | Bitnami Neo4j AMI ID.                            |
| APOC_VERSION                 | 4.4.0.3                            | string | APOC version for Neo4j.                          |
| GDS_VERSION                  | 2.0.1                              | string | GDS version for Neo4j.                           |
| GITHUB_PERSONAL_ACCESS_TOKEN | <secret value>                     | string | GitHub PAT for repository access.                |

***Important**:* *Always use a `.env` file or AWS SSM Parameter Store or Secrets Manager for sensitive variables like credentials and API keys. Never hard-code them, including when developing. AWS will quarantine an account if any credentials get accidentally exposed and this will cause problems. Make sure to update `.gitignore` to avoid pushing sensitive data to public repositories.*

### Makefile Usage
Once an AWS profile is configured and environment variables are exported, the application can be deployed using `make`. You are required to specify the `STAGE` variable everytime `make` is called to ensure that the correct environment is selected when there are multiple deployments.
```bash
STAGE=<stage> make deploy
```
It is also possible to deploy or update the database or pipeline services.
```bash
# Deploy/update only the database service
STAGE=<stage> make database.deploy

# Deploy/update only the pipeline service
STAGE=<stage> make pipeline.deploy

# Deploy/update only the pipeline serverless stack
STAGE=<stage> make pipeline.functions.deploy

# Deploy/update only the Docker image for the build job
STAGE=<stage> make pipeline.jobs.deploy
```
*Note:* It is recommended to only deploy from the project root. This is because common parameters are passed from the root Makefile to nested Makefiles. If a stack has not been changed, the deployment script will continue until it reaches a stack with changes and deploy that.

#### Makefile Command Reference
To see a list of possible commands using Make, run `make` on the command line. You can also refer to the `Makefile Usage` section in the [Sphinx documentation](#documentation).
```bash
# Deploy all CloudFormation based services
STAGE=<stage> make deploy

# Deploy config files and scripts to S3
STAGE=<stage> make config.deploy

# Run the StepFunctions State Machine to load Neo4j
STAGE=<stage> make database.load.run releases=<version> align=<boolean> kir=<boolean> limit=<int>

# Retrieve Neo4j credentials after deployment
STAGE=<stage> make database.get.credentials

# Retrieve Neo4j URL after deployment
STAGE=<stage> make database.get.endpoint

# Download logs from EC2
STAGE=<stage> make get.logs

# Download CSV data from S3
STAGE=<stage> make get.data

# Delete all CloudFormation based services and data, default is data=false
STAGE=<stage> make delete data=<true/false>

# Delete a specific layer
STAGE=<stage> make pipeline.delete

# Subscribe an email for notifications (unsubscribe using console)
STAGE=<stage> make monitoring.subscribe-email email=<email>
```

## Managing Configuration
Configuration is managed using JSON files, SSM Parameter Store, Secrets Manager, and shell variables. To deploy changes in these files, run the command.
```bash
STAGE=<stage> make config.deploy
```

### Database Configuration

#### Neo4j
Custom configuration settings for Neo4j are contained in `neo4j.template`. This file is copied into `/etc/neo4j` during boot or can be done manually. When Neo4j is restarted it will use the settings in `neo4j.template` to overwrite `neo4j.conf`. More information can be found in the documentation here: [Neo4j Cloud Virtual Machines](https://neo4j.com/developer/neo4j-cloud-vms/).

#### Shell Scripts
Bash scripts are used for automating Neo4j configuration, loading and backup. These are stored in S3 and run using SSM Run Command. These are found in `gfe-db/gfe-db/database/scripts/`.

```bash
gfe-db/database/scripts
├── Makefile                  # Orchestrates tasks on the database instance
├── init  
│   ├── create_cert.sh        # Create an SSL certificate
│   └── eip_assoc_waiter.sh   # Waits for the instance to associate with an Elastic IP
├── load_db.sh                # Loads data into Neo4j
├── send_heartbeat.sh         # Sends task heartbeat to Step Functions API during loading
└── start_task.sh             # Coordinates database loading with the Step Functions API
```

To update shell scripts on the Neo4j instance, run the following command.
```bash
# sync the scripts from S3 to the instance (using Systems Manager Run Command)
STAGE=<stage> make database.sync-scripts
```

#### Cypher Scripts
Cypher scripts manage node constraints & indexes and load the data. These are found in `gfe-db/gfe-db/database/neo4j/cypher/`.

```bash
gfe-db/database/neo4j/
├── cypher                      # Cypher scripts
│   ├── create_constraints.cyp  # Creates constraints and indexes
│   ├── drop_constraints.cyp    # Drops constraints and indexes
│   ├── init.cyp                # Run intitialization queries
│   └── load.cyp                # Load Neo4j from local files
├── (neo4j.conf)                # Artifact create by Makefile containing DNS configuration
└── neo4j.template              # Config template used by Makefile
```

#### SSL Policy
SSL is created during deployment and requires access on ports 443 (HTTPS) and 7473 (Neo4j browser for HTTPS).
### Data Pipeline Config

#### Invocation Input
Base input parameters (excluding the `releases` value) are passed to the Step Functions State Machine and determine it's behavior during build. The `releases` value is appended at runtime by the trigger Lambda when it finds a new release in the source repository. The `pipeline-input.json` is stored in S3 and contains the default configuration used for automated updates.
```json
// pipeline-input.json
{
  "align": "False",
  "kir": "False",
  "mem_profile": "False",
  "limit": ""
}

```
| Variable    | Example Value | Type   | Description                                                        |
| ----------- | ------------- | ------ | ------------------------------------------------------------------ |
| LIMIT       | 1000          | string | Number of alleles to build. Leave blank ("") to build all alleles. |
| ALIGN       | False         | string | Include or exclude alignments in the build                         |
| KIR         | False         | string | Include or exclude KIR data alignments in the build                |
| MEM_PROFILE | False         | string | Enable memory profiling (for catching memory leaks during build)   |

The data pipeline can also be invoked from the command line:
```bash
STAGE=<stage> make database.load.run releases=<version> align=<boolean> kir=<boolean> limit=<int>
```

#### IMGT/HLA Release Versions State
The application's state tracks which releases have been processed and added to the database. This file tracks the releases which have already been processed. If the `gfe-db-invoke-pipeline` function detects a valid release branch in the source data repository that is not in the `releases` array, it will start the pipeline for this release. Once the update is finished, the processed release is appended to the array.
```json
// ./gfe-db/gfe-db/pipeline/config/IMGTHLA-repository-state.json
{
  "timestamp": "2021-12-09 02:36:59",
  "repository_url": "https://github.com/ANHIG/IMGTHLA",
  // Releases already loaded (or skipped)
  // Any releases not found in this array will be loaded
  "releases": [
    "3100",
    ...,
    "3510"
  ]
}
```

| Variable       | Example Value                    | Type             | Description                                                                                                               |
| -------------- | -------------------------------- | ---------------- | ------------------------------------------------------------------------------------------------------------------------- |
| repository_url | https://github.com/ANHIG/IMGTHLA | string           | The repository the trigger is watching                                                                                    |
| releases       | ["3100", ..., "3510"]            | array of strings | List of available releases. Any release added to the repository that is not in this list will trigger the pipeline build. |

## Loading Neo4j
For each invocation the data pipeline will download raw data from [ANHIG/IMGTHLA](https://github.com/ANHIG/IMGTHLA) GitHub repository, build a set of intermediate CSV files and load these into Neo4j via S3. To invoke the pipeline, run the following command.
```bash
STAGE=<stage> make database.load.run releases="<version>"

# Example for single version
STAGE=<stage> make database.load.run releases="3510"

# Example for multiple versions
STAGE=<stage> make database.load.run releases="3490,3500,3510"

# Example with limit
STAGE=<stage> make database.load.run releases="3510" limit="1000"

# Example with all arguments included
STAGE=<stage> make database.load.run releases="3510" limit="" align="False" kir="False"
```

These commands build an event payload to send to the `invoke-gfe-db-pipeline` Lambda.
```json
// Test payload example
{
  "align": false,
  "kir": false,
  "mem_profile": false,
  "limit": "",
  "releases": 3510
}
```

The Lambda function returns the following object which can be viewed in CloudWatch Logs.
```json
// invoke-gfe-db-pipeline Lambda function response
{
  "status": 200,
  "message": "Pipeline triggered",
  "input": [
    {
      "ALIGN": "False",
      "KIR": "False",
      "MEM_PROFILE": "False",
      "LIMIT": "",
      "RELEASES": "3510"
    },
    ...
  ]
}
```

### Clean Up
To tear down resources run the command. You will need to manually delete the data in the S3 bucket first to avoid an error in CloudFormation.
```bash
STAGE=<stage> make delete data=<true/false>
```
Use the following commands to tear down individual services. Make sure to [backup](#backup--restore) your data first.
```bash
# Delete only the database service
STAGE=<stage> make database.delete

# Delete only the pipeline service
STAGE=<stage> make pipeline.delete
```

## Backup & Restore

### Backups

Backups are orchestrated by Systems Manager and automated everyday at midnight US/Central time by default. To create a backup, run the command.

```bash
STAGE=<stage> make database.backup
```

This will create a backup of the Neo4j database and store it in S3 under the path `s3://<data bucket name>/backups/neo4j/YYYY/MM/DD/HH/gfedb.zip`.

### Restore

To see a list of available backup dates that can be restored, run the command.

```bash
STAGE=<stage> make database.backup.list
```

To restore from a backup, pass the date of the backup you wish to restore using the format YYYY/MM/DD/HH.

```bash
STAGE=<stage> make database.restore from_date=<YYYY/MM/DD/HH>
```

## Local Development

### Creating a Python Virtual Environment
When developing locally, you will need to create an individual virtual environment to run scripts in the `jobs` or `functions` directories, since they require different dependencies.
```bash
python3 -m venv .venv-dev
source .venv-dev/bin/activate
pip install -U pip
pip install -r requirements.txt
```

To use the virtual environment inside a Jupyter Notebook, first activate the virtual environment, then create a kernel for it.
```bash
# Install ipykernal
pip install ipykernel python-dotenv

# Add the kernel
python3 -m ipykernel install --user --name=<environment name>

# Remove the kernel
jupyter kernelspec uninstall <environment name>
```

## Documentation
It is not necessary to install Sphinx to view `gfe-db` documentation because it is already built and available in the `docs/` folder, but you will need it to edit them. To get the local `index.html` path run the command and navigate to the URL in a browser.

```bash
STAGE=<stage> make docs.url
```

### Editing and Building the Documentation
Create a virtual environment and install the dependencies for working with the documentation.
```bash
python3 -m venv .venv-docs
source .venv-docs/bin/activate
pip install -U pip
pip install -r requirements-docs.txt
```

After making your edits, you can build the HTML assets by running the command.
```bash
STAGE=<stage> make docs.build
```

## Troubleshooting
* Check your AWS credentials in `~/.aws/credentials`
* Check that the correct environment variables are set in `.env`
* Check that Python 3.9 is being used
* Make sure you are accessing Neo4j Browser using HTTPS, not HTTP (some browsers like Chrome will not show `https://` in the URL making it hard to tell)

## Authors
**Primary Contact:** Martin Maiers ([@mmaiers-nmdp](https://github.com/mmaiers-nmdp))\
**Contact:** Pradeep Bashyal ([@pbashyal-nmdp](https://github.com/pbashyal-nmdp))\
**Contact:** Gregory (Chris) Lindsey ([@chrisammon3000](https://github.com/chrisammon3000))

<!-- TODO make sure these are up to date -->
## References & Links
 * [bioinformatics.bethematchclinical.org](https://bioinformatics.bethematchclinical.org)
 * [Cypher Reference](https://neo4j.com/docs/cypher-manual/current/)
 * [Getting Certificates for Neo4j with LetsEncrypt](https://medium.com/neo4j/getting-certificates-for-neo4j-with-letsencrypt-a8d05c415bbd)

-----------------
<br>
<p align="center">
  <img src="https://bethematch.org/content/site/images/btm_logo.png" alt="Be The Match>
</p>
