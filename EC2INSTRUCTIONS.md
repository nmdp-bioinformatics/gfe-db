# EC2 Deployment

This covers how to build and deploy `gfe-db` on an EC2 instance running Neo4j in Docker.

## Table of Contents
<!-- - [EC2 Deployment](#ec2-deployment)
  - [Table of Contents](#table-of-contents) -->
  - [Setup EC2](#setup-ec2)
    - [Security Groups](#security-groups)
      <!-- - [Neo4j Access](#neo4j-access)
      - [SSH Access](#ssh-access) -->
    - [Launch Instance](#launch-instance)
    - [SSH into EC2](#ssh-into-ec2)
  - [Installation](#installation)
    - [Dependencies](#dependencies)
    - [Docker](#docker)
    - [gfe-db](#gfe-db)
    - [Create Virtual Environment](#create-virtual-environment)
    - [Export Environment Variables](#export-environment-variables)
  - [Run Neo4j](#run-neo4j)
    - [Configure Neo4j for EC2](#configure-neo4j-for-ec2)
    - [Build Neo4j Docker](#build-neo4j-docker)
    - [Start Neo4j](#start-neo4j)
  - [Build gfe-db](#build-gfe-db)
- [References & Links](#references--links)

## Setup EC2

### Security Groups
#### Neo4j Access
In the AWS EC2 Network & Security console or with AWS CLI create a security group for Neo4j and name it `neo4j-sg` or something similar. Leave the outbound rules as their default values and configure the inbound rules according to the table.
| Type       | Protocol | Port Range | Source    | Description            |
|------------|----------|------------|-----------|------------------------|
| Custom TCP | TCP      | 7474       | 0.0.0.0/0 | HTTP for Neo4j         |
| Custom TCP | TCP      | 7473       | 0.0.0.0/0 | HTTPS access for Neo4j |
| Custom TCP | TCP      | 7687       | 0.0.0.0/0 | Bolt access for Neo4j  |

#### SSH Access
Create another security group for SSH access and name it `gfe-db-ssh` or something similar. Leave the outbound rules as their default values and configure the inbound rules according to the table. Select "My IP" for the source.
| Type | Protocol | Port Range | Source    | Description                |
|------|----------|------------|-----------|----------------------------|
| SSH  | TCP      | 22         | \<your IP> | SSH access to Neo4j server|

### Launch Instance
* Using the AWS EC2 console or AWS CLI, launch an instance with the following configuration.
* Amazon Linux 2 AMI
* `c5d.18xlarge`
* Enable both **Auto-assign Public IP** and **Enable termination protection**.
* No EBS volume is needed for instances with "d" in the name (`c5d... etc.`)
* Give it a tag with key "Name" and value "gfe-db"
* Select the `neo4j-sg` and `gfe-db-ssh` security groups
* Create a key or select an existing one

### SSH into EC2
Make sure the key you are using has the correct permissions before using it for SSH access.
```bash
chmod 400 <your-key.pem>
```

When the instance is ready, copy the Public IPv4 DNS and ssh into it with this command.
```bash
# SSH into the instance
ssh -i <your-key.pem> ec2-user@<Public IPv4 DNS>
```

## Installation
### Dependencies
This step installs git and Python 3.8 on EC2.
```bash
# Update the instance
sudo yum update -y

# Install git
sudo yum install git -y

# Install python 3.8
sudo amazon-linux-extras install python3.8
```

### Docker
This step installs Docker with the following commands.
```bash
# Install Docker
sudo amazon-linux-extras install docker

# Start the Docker service
sudo service docker start

# Add the current ec2-user to the docker group 
sudo usermod -a -G docker ec2-user
```
Log out and log back in to the instance.
```bash
# Exit the shell
exit

# SSH into the instance
ssh -i <key_pair_name> ec2-user@<Public_IPv4_DNS>
```
Confirm that Docker is running. It should list the client and server details.
```bash
docker info
```

### gfe-db
Clone the repo.
```bash
git clone https://github.com/abk7777/gfe-db.git
```
### Create Virtual Environment
Create a virtual environment and activate. Make sure you are using Python 3.8.
```bash
# Create .venv specifying Python 3.8
python3.8 -m venv .venv

# Activate
source .venv/bin/activate
```
Install the requirements.
```bash
pip install -r requirements.txt
```
### Export Environment Variables
Make sure the environment variables in `bin/set_env.sh ` are exported to the environment.
```bash
source bin/set_env.sh 
```

## Run Neo4j
### Configure Neo4j for EC2
In the Dockerfile, update the following settings using nano or vim.
```docker
ENV NEO4J_dbms_memory_heap_initial__size=32G
ENV NEO4J_dbms_memory_heap_max__size=32G
ENV NEO4J_dbms_memory_pagecache_size=4G
```
### Build Neo4j Docker
Build the Docker image as defined in the Dockerfile. 
```
docker build --tag gfe-db .
```
### Start Neo4j
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
Neo4j will change permissions on the `data/` directory. You will need to update them before the build step.
```bash
# Update permissions from project root
cd gfe-db
sudo -R chmod 777 .
```

## Build gfe-db
Run this script to generate a set CSV files of GFE data in the `data/csv/` directory. It is recommended to limit the number of alleles and start with a small number to avoid excessive build and load times.
```bash
# Limit the build to 1000 alleles (recommended for testing)
bash bin/build.sh 1000

# Build complete database (takes a while)
bash bin/build.sh
```



# References & Links
* [AWS Quick Start Guide: Launch a Linux Virtual Machine](https://docs.aws.amazon.com/quickstarts/latest/vmlaunch/welcome.html)
* [Docker basics for Amazon ECS](https://docs.aws.amazon.com/AmazonECS/latest/developerguide/docker-basics.html) (Docker installation steps apply to EC2)
* [Neo4j Performance Tuning](https://neo4j.com/developer/guide-performance-tuning/)

<br>
<p align="center">
  <img src="https://bethematch.org/content/site/images/btm_logo.png">
</p>