# This Makefile only orchestrates process that are run on the EC2 database instance, it is not called by the root Makefile
# Run as sudo

# Application specific environment variables
include .env
export

export ROOT_DIR=/home/bitnami
export NEO4J_HOME=/home/bitnami/neo4j

target:
	$(info ${HELP_MESSAGE})
	@exit 0

# TODO fetch all shell scripts from S3
# fetch.scripts:
# fetch.cypher

# env.fetch:

env.check:
ifndef DATA_BUCKET_NAME
$(error DATA_BUCKET_NAME is not set. Please select an AWS profile to use.)
endif
ifndef NEO4J_HOME
$(error NEO4J_HOME is not set.)
endif
ifndef HOST_DOMAIN
$(error HOST_DOMAIN is not set.)
endif
ifndef ADMIN_EMAIL
$(error ADMIN_EMAIL is not set.)
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Found environment variables" 2>&1

env.app.list:
	@cat .env

env.make.list:
	@printenv

# # TODO add dependency checks
# deps.check:
# 	# awscli, s3 connection
# 	# jq
# 	# curl
# 	# wget
# 	# certbot
# 	# python

# aws.cfn.helpers.install

# aws.cfn.cfn-signal.success

# aws.cfn.cfn-signal.failure

# aws.ssm.agent.install

# aws.cw.agent.install

# eip
eip.waiter:
	@bash ${ROOT_DIR}/user_data/eip_assoc_waiter.sh

# neo4j
neo4j.stop:
	@/opt/bitnami/ctlscript.sh stop

neo4j.start:
	@/opt/bitnami/ctlscript.sh start

neo4j.restart:
	@/opt/bitnami/ctlscript.sh restart

# TODO test
neo4j.config.update:
	@mv ${NEO4J_HOME}/conf/neo4j.conf ${NEO4J_HOME}/conf/neo4j.conf.bkp
	@echo "Downloading Neo4j configuration from ${DATA_BUCKET_NAME}"
	@aws s3 cp s3://${DATA_BUCKET_NAME}/config/neo4j/neo4j.conf ${NEO4J_HOME}/conf/neo4j.conf

# neo4j.plugins.update.apoc:

# TODO test
# TODO add version argument
neo4j.plugins.update.gds:
	@echo "Downloading Neo4j Graph Data Science plugin..."
	@curl -L https://graphdatascience.ninja/neo4j-graph-data-science-2.0.1.zip -O
	@unzip neo4j-graph-data-science-2.0.1.zip
	@mv neo4j-graph-data-science-2.0.1.jar $NEO4J_HOME/plugins/neo4j-graph-data-science-2.0.1.jar

# TODO test
# neo4j.graph.init:
# # Run cypher init scripts
# mkdir -p $NEO4J_HOME/{backups,cypher}
# echo "Fetching cypher scripts from S3..."
# aws s3 cp --recursive s3://${DataBucketName}/config/neo4j/cypher/ $NEO4J_HOME/cypher/
# cat $NEO4J_HOME/cypher/${Neo4jInitScript} | \
# /usr/bin/cypher-shell \
#   -u neo4j \
#   -p ${Neo4jPassword}

# TODO test
# neo4j.graph.import:
# mkdir -p $NEO4J_HOME/backups $NEO4J_HOME/cypher
# echo "Fetching cypher scripts from S3..."
# aws s3 cp --recursive s3://${DataBucketName}/config/neo4j/cypher/ $NEO4J_HOME/cypher/
# cat $NEO4J_HOME/cypher/${Neo4jInitScript} | \
# /usr/bin/cypher-shell \
#   -u neo4j \
#   -p ${Neo4jPassword}

# neo4j.graph.constraints.create:

# neo4j.graph.constraints.delete:

# TODO test
# ssl.create
	@bash user_data/create_cert.sh ${HostDomain} ${AdminEmail}

# ssl.renew

# TODO test
# logs
# echo "Copying logs to S3..."
# export INSTANCE_ID=$(curl -H "X-aws-ec2-metadata-token: $TOKEN" -v http://169.254.169.254/latest/meta-data/instance-id)
# journalctl -b > /tmp/logs/$INSTANCE_ID-boot.log
# echo "Copying Neo4j logs to S3..."
# journalctl -e -u neo4j > /tmp/logs/neo4j.log
# aws s3 cp --recursive /tmp/logs/ s3://${DataBucketName}/logs/database/bootstrap/$(date +'%Y/%m/%d/%H')/
# echo "Ready"

define HELP_MESSAGE

	Environment variables:

	SERVICE: "${SERVICE}"
		Description: Name of the service being deployed

	Common usage:

	...::: Deploy all CloudFormation based services :::...
	$ make deploy

	...::: Delete all CloudFormation based services :::...
	$ make delete

endef