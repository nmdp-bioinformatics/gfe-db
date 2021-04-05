# EC2 Deployment

This covers how to build and deploy `gfe-db` on an EC2 instance running Neo4j in Docker.

## Setup EC2

### Create Security Groups
In the AWS EC2 Network & Security console or with AWS CLI create a security group for Neo4j and name it `neo4j-sg` or something similar. Leave the outbound rules as their default values and configure the inbound rules according to the table.
| Type       | Protocol | Port Range | Source    | Description            |
|------------|----------|------------|-----------|------------------------|
| Custom TCP | TCP      | 7474       | 0.0.0.0/0 | HTTP for Neo4j         |
| Custom TCP | TCP      | 7473       | 0.0.0.0/0 | HTTPS access for Neo4j |
| Custom TCP | TCP      | 7687       | 0.0.0.0/0 | Bolt access for Neo4j  |

Create another security group for SSH access and name it `gfe-db-ssh` or something similar. Leave the outbound rules as their default values and configure the inbound rules according to the table. Select "My IP" for the source.
| Type | Protocol | Port Range | Source    | Description                |
|------|----------|------------|-----------|----------------------------|
| SSH  | TCP      | 22         | \<your IP> | SSH access to Neo4j server|

### Launch Instance
* Using the AWS EC2 console or AWS CLI, launch an instance with the following configuration.
* Amazon Linux 2 AMI
* `c5d.18xlarge`
* Enable both **Auto-assign Public IP** and **Enable termination protection**.
* No EBS volume (NVMe SSD will be used instead for any instance ending )
* Give it a tag with key `Name` and value `gfe-db`
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

## Install Docker
Inside the instance install Docker with the following commands.
```bash
# Update the instance
sudo yum update -y

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

## Configure Neo4j for EC2
In the Dockerfile, update the following settings.
```docker
ENV NEO4J_dbms_memory_heap_initial__size=32G
ENV NEO4J_dbms_memory_heap_max__size=32G
ENV NEO4J_dbms_memory_pagecache_size=4G
```
Continue with the build and load process as described in the main README.

# References & Links
* [AWS Quick Start Guide: Launch a Linux Virtual Machine](https://docs.aws.amazon.com/quickstarts/latest/vmlaunch/welcome.html)
* [Docker basics for Amazon ECS](https://docs.aws.amazon.com/AmazonECS/latest/developerguide/docker-basics.html) (Docker installation steps apply to EC2)
* [Neo4j Performance Tuning](https://neo4j.com/developer/guide-performance-tuning/)

<br>
<p align="center">
  <img src="https://bethematch.org/content/site/images/btm_logo.png">
</p>