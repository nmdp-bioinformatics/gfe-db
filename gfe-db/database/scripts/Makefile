# This Makefile only orchestrates process that are run on the EC2 database instance, it is deployed and called on the server
# Run as sudo only

# Application specific environment variables
include env.sh
export

export BITNAMI_HOME ?= /home/bitnami
export LOGS_DIR := /tmp/logs
export NEO4J_HOME := /opt/bitnami/neo4j
export NEO4J_USERNAME := $(shell cat ${BITNAMI_HOME}/bitnami_credentials | sed 's/ /\n/g' | grep -E "'([a-zA-Z0-9]+)'" | sed "s/'//g" | sed "s/\.//g" | head -n 1)
export NEO4J_PASSWORD := $(shell cat ${BITNAMI_HOME}/bitnami_credentials | sed 's/ /\n/g' | grep -E "'([a-zA-Z0-9]+)'" | sed "s/'//g" | sed "s/\.//g" | tail -n 1)
export INSTANCE_ID := $(shell curl -sS http://169.254.169.254/latest/meta-data/instance-id)

target:
	$(info ${HELP_MESSAGE})
	@exit 0

env:
	@printenv

# TODO fetch all shell scripts from S3
# get.scripts:
# get.cypher:
get.data: #=> Download the build data locally; get.data release=3470
	# should be [[ -n "$$release" ]]
	@[ -z "$$release" ] && aws s3 cp --recursive s3://${DATA_BUCKET_NAME}/data/$$release/csv/ ${NEO4J_HOME}/import/ || echo "No release argument"
# env.fetch:

env.check:
ifndef DATA_BUCKET_NAME
$(error DATA_BUCKET_NAME is not set.)
endif
ifndef NEO4J_HOME
$(error NEO4J_HOME is not set.)
endif
ifndef NEO4J_USERNAME
$(error NEO4J_HOME is not set.)
endif
ifndef NEO4J_PASSWORD
$(error NEO4J_HOME is not set.)
endif
ifndef HOST_DOMAIN
$(error HOST_DOMAIN is not set.)
endif
ifndef ADMIN_EMAIL
$(error ADMIN_EMAIL is not set.)
endif
ifndef APOC_VERSION
$(error APOC_VERSION is not set.)
endif
ifndef GDS_VERSION
$(error GDS_VERSION is not set.)
endif
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Found environment variables" 2>&1

# # TODO add dependency checks
# deps.check:
# 	# awscli, s3 connection
# 	# jq
# 	# curl
# 	# wget
# 	# certbot
# 	# python

bootstrap:
	# $(MAKE) cfn-helpers.install
	# $(MAKE) ssm.install
	# $(MAKE) cloudwatch.install
	# $(MAKE) cfn-signal exitcode=0
	$(MAKE) eip.waiter
	$(MAKE) ssl.create-cert

neo4j: bootstrap
	$(MAKE) neo4j.update-credentials
	$(MAKE) neo4j.config.update
	$(MAKE) neo4j.plugins.install-apoc
	$(MAKE) neo4j.plugins.install-gds
	$(MAKE) neo4j.restart
	$(MAKE) neo4j.init
	$(MAKE) neo4j.restart
	$(MAKE) copy-logs

# TODO keep in User Data to separate the boostrap logic from the application logic 
# cfn-helpers.install:
# 	@mkdir -p /tmp/cfn
# 	@wget -q -O /tmp/cfn/aws-cfn-bootstrap-py3-latest.tar.gz https://s3.amazonaws.com/cloudformation-examples/aws-cfn-bootstrap-py3-latest.tar.gz
# 	@python3 -m easy_install --script-dir /opt/aws/bin /tmp/cfn/aws-cfn-bootstrap-py3-latest.tar.gz

# TODO keep in User Data to separate the boostrap logic from the application logic 
# cfn-signal: ##=> make cfn-signal exitcode=0 for success
# 	@res=$$(/opt/aws/bin/cfn-signal \
# 		--exit-code $$exitcode \
# 		--stack  ${STAGE}-${APP_NAME}-database \
# 		--resource Neo4jDatabaseInstance \
# 		--region ${AWS_REGION}) || echo "skipping cfn-signal"

# TODO keep in User Data to separate the boostrap logic from the application logic 
# ssm.install:
# 	@mkdir -p /tmp/ssm
# 	@wget -q -O /tmp/ssm/amazon-ssm-agent.deb https://s3.amazonaws.com/ec2-downloads-windows/SSMAgent/latest/debian_amd64/amazon-ssm-agent.deb
# 	@dpkg -i /tmp/ssm/amazon-ssm-agent.deb
# 	@systemctl enable amazon-ssm-agent

# TODO keep in User Data to separate the boostrap logic from the application logic 
# cloudwatch.install:
# 	@mkdir -p /tmp/amazon-cloudwatch-agent
# 	@wget -q -O /tmp/amazon-cloudwatch-agent/amazon-cloudwatch-agent.deb https://s3.amazonaws.com/amazoncloudwatch-agent/debian/amd64/latest/amazon-cloudwatch-agent.deb
# 	@dpkg -i -E /tmp/amazon-cloudwatch-agent/amazon-cloudwatch-agent.deb
# 	@aws s3 cp s3://${DATA_BUCKET_NAME}/config/amazon-cloudwatch-agent/amazon-cloudwatch-agent.json /opt/aws/amazon-cloudwatch-agent/etc/amazon-cloudwatch-agent.json
# 	@systemctl enable amazon-cloudwatch-agent.service
# 	@service amazon-cloudwatch-agent start

# eip
eip.waiter:
	@bash ${BITNAMI_HOME}/init/eip_assoc_waiter.sh

# neo4j
neo4j.update-credentials:
	@aws --region ${AWS_REGION} secretsmanager update-secret \
		--secret-id ${NEO4J_CREDENTIALS_SECRET_ARN} \
		--secret-string "{\"NEO4J_USERNAME\":\"${NEO4J_USERNAME}\",\"NEO4J_PASSWORD\":\"${NEO4J_PASSWORD}\"}"

# TODO update database env variables using SSM
# neo4j.update-env:
# 	@echo "STAGE=${Stage}" >> ${BITNAMI_HOME}/env.sh
# 	@echo "APP_NAME=${AppName}" >> ${BITNAMI_HOME}/env.sh
# 	@echo "AWS_REGION=${AWS::Region}" >> ${BITNAMI_HOME}/env.sh
# 	@echo "DATA_BUCKET_NAME=${DataBucketName}" >> ${BITNAMI_HOME}/env.sh
# 	@echo "HOST_DOMAIN=${HostDomain}" >> ${BITNAMI_HOME}/env.sh
# 	@echo "SUBDOMAIN=${Subdomain}" >> ${BITNAMI_HOME}/env.sh
# 	@echo "ADMIN_EMAIL=${AdminEmail}" >> ${BITNAMI_HOME}/env.sh
# 	@echo "APOC_VERSION=${APOCVersion}" >> ${BITNAMI_HOME}/env.sh
# 	@echo "GDS_VERSION=${GDSVersion}" >> ${BITNAMI_HOME}/env.sh
# 	@echo "BITNAMI_HOME=/home/bitnami" >> ${BITNAMI_HOME}/env.sh
# 	@echo "NEO4J_HOME=/opt/bitnami/neo4j" >> ${BITNAMI_HOME}/env.sh
# 	@echo "BITNAMI_NEO4J=/bitnami/neo4j" >> ${BITNAMI_HOME}/env.sh
# 	@echo "NEO4J_CREDENTIALS_SECRET_ARN=${Neo4jCredentialsSecret}" >> ${BITNAMI_HOME}/env.sh
# 	@echo "# Make application variables available" >> ${BITNAMI_HOME}/.bashrc
# 	@echo "set -a && source /home/bitnami/env.sh && set +a" >> ${BITNAMI_HOME}/.bashrc

neo4j.status:
	@/opt/bitnami/ctlscript.sh status

neo4j.stop:
	@/opt/bitnami/ctlscript.sh stop

neo4j.start:
	@/opt/bitnami/ctlscript.sh start
	$(MAKE) neo4j.waiter

neo4j.restart:
	@/opt/bitnami/ctlscript.sh restart
	$(MAKE) neo4j.waiter

neo4j.config.update:
	@mv ${NEO4J_HOME}/conf/neo4j.conf ${NEO4J_HOME}/conf/neo4j.conf.bkp
	@echo "Downloading Neo4j configuration from ${DATA_BUCKET_NAME}"
	@aws s3 cp s3://${DATA_BUCKET_NAME}/config/neo4j/neo4j.conf ${NEO4J_HOME}/conf/neo4j.conf

neo4j.plugins.install-apoc:
	@echo "Downloading APOC Full plugin..."
	@curl -L https://github.com/neo4j-contrib/neo4j-apoc-procedures/releases/download/${APOC_VERSION}/apoc-${APOC_VERSION}-all.jar -O
	@mv apoc-${APOC_VERSION}-all.jar ${NEO4J_HOME}/plugins/apoc-${APOC_VERSION}-all.jar

neo4j.plugins.install-gds:
	@echo "Downloading Neo4j Graph Data Science plugin..."
	@curl -L https://graphdatascience.ninja/neo4j-graph-data-science-${GDS_VERSION}.zip -O
	@unzip neo4j-graph-data-science-${GDS_VERSION}.zip
	@mv neo4j-graph-data-science-${GDS_VERSION}.jar ${NEO4J_HOME}/plugins/neo4j-graph-data-science-${GDS_VERSION}.jar
	@echo "Cleaning up"
	@rm neo4j-graph-data-science-${GDS_VERSION}.zip

neo4j.init:
	@mkdir -p ${NEO4J_HOME}/backups ${NEO4J_HOME}/cypher
	@echo "Fetching Cypher scripts from S3..."
	@aws s3 cp --recursive s3://${DATA_BUCKET_NAME}/config/neo4j/cypher/ ${NEO4J_HOME}/cypher/
	@echo "Executing initialization queries..."
	@cat ${NEO4J_HOME}/cypher/init.cyp | /opt/bitnami/neo4j/bin/cypher-shell -u ${NEO4J_USERNAME} -p ${NEO4J_PASSWORD} -a neo4j+s://${SUBDOMAIN}.${HOST_DOMAIN}:7687

neo4j.waiter:
	@timeout=120 && \
	counter=0 && \
	echo "Waiting for response from Neo4j at https://${SUBDOMAIN}.${HOST_DOMAIN}:7473..." && \
	until $$(curl --output /dev/null --silent --head --fail https://${SUBDOMAIN}.${HOST_DOMAIN}:7473) ; do \
		printf '.' ; \
		sleep 1 ; \
		counter=$$((counter + 1)) ; \
		[ $$counter -eq $$timeout ] && break || true ; \
	done && \
	printf "%s\n" " " && \
	[ $$counter -eq $$timeout ] && echo "Operation timed out!" || echo "Neo4j is ready"

# TODO add target for loading
# neo4j.load: # params=<json string>
# 	@bash start_task.sh $$params

neo4j.backup:
	$(MAKE) neo4j.stop
	@cd ${BITNAMI_NEO4J}/ && zip -r ${BITNAMI_HOME}/gfedb.zip data
	@result=$$(aws s3 cp ${BITNAMI_HOME}/gfedb.zip s3://${DATA_BUCKET_NAME}/backups/neo4j/$$(date +'%Y/%m/%d/%H')/gfedb.zip) && \
		echo $$result
	@rm ${BITNAMI_HOME}/gfedb.zip
	$(MAKE) neo4j.start

# TODO add server-side validation before running restore
	# @echo "test: $$from_date"
	# @[ "$$from_date" != "" ] || (echo "from_date is required" && exit 1) && \
	# [ "$$(echo "$$from_date" | grep -E '^[0-9]{4}/[0-9]{2}/[0-9]{2}/[0-9]{2}$$')" ] || (echo "from_date must be in the format YYYY/MM/DD/HH" && exit 1) && \
	# echo "running restore for $$from_date"
neo4j.restore: #from_date=<YYYY/MM/DD/HH>
	$(MAKE) neo4j.stop
	@aws s3 cp s3://${DATA_BUCKET_NAME}/backups/neo4j/$$from_date/gfedb.zip ${BITNAMI_HOME}/gfedb.zip
	@unzip -o ${BITNAMI_HOME}/gfedb.zip
	@cp -r data/* ${BITNAMI_NEO4J}/data/
	@chown -R neo4j:neo4j ${BITNAMI_NEO4J}/data/
	@rm ${BITNAMI_HOME}/gfedb.zip
	# @rm -r ${BITNAMI_HOME}/data/
	$(MAKE) neo4j.start

# # TODO test
# neo4j.import:
# 	@mkdir -p ${NEO4J_HOME}/backups ${NEO4J_HOME}/cypher
# 	@echo "Fetching cypher scripts from S3..."
# 	@aws s3 cp --recursive s3://${DATA_BUCKET_NAME}/config/neo4j/cypher/ ${NEO4J_HOME}/cypher/
# 	@cat ${NEO4J_HOME}/cypher/${Neo4jInitScript} | cypher-shell -u neo4j -p ${NEO4J_PASSWORD}

neo4j.create-constraints:
	@mkdir -p ${NEO4J_HOME}/backups ${NEO4J_HOME}/cypher
	@echo "Fetching cypher scripts from S3..."
	@aws s3 cp --recursive s3://${DATA_BUCKET_NAME}/config/neo4j/cypher/ ${NEO4J_HOME}/cypher/
	@cat ${NEO4J_HOME}/cypher/create_constraints.cyp | cypher-shell -u ${NEO4J_USERNAME} -p ${NEO4J_PASSWORD} -a neo4j+s://${SUBDOMAIN}.${HOST_DOMAIN}:7687

neo4j.drop-constraints:
	@mkdir -p ${NEO4J_HOME}/backups ${NEO4J_HOME}/cypher
	@echo "Fetching cypher scripts from S3..."
	@aws s3 cp --recursive s3://${DATA_BUCKET_NAME}/config/neo4j/cypher/ ${NEO4J_HOME}/cypher/
	@cat ${NEO4J_HOME}/cypher/drop_constraints.cyp | cypher-shell -u ${NEO4J_USERNAME} -p ${NEO4J_PASSWORD} -a neo4j+s://${SUBDOMAIN}.${HOST_DOMAIN}:7687

ssl.create-cert:
	@bash init/create_cert.sh "${SUBDOMAIN}.${HOST_DOMAIN}" ${ADMIN_EMAIL}

copy-logs:
	@mkdir -p ${LOGS_DIR}/neo4j ${LOGS_DIR}/system
	$(MAKE) logs.bootstrap
	$(MAKE) logs.system
	$(MAKE) logs.neo4j
	@echo "Copying logs to S3..." && \
	aws s3 cp --recursive /tmp/logs/ s3://${DATA_BUCKET_NAME}/logs/database/${INSTANCE_ID}/$$(date +'%Y/%m/%d/%H')/

logs.bootstrap:
	@journalctl -b > ${LOGS_DIR}/bootstrap.log

logs.system:
	@cp -r /var/log/* ${LOGS_DIR}/system

logs.neo4j:
	@cp -r ${NEO4J_HOME}/logs/* ${LOGS_DIR}/neo4j/

define HELP_MESSAGE

	Environment variables:

	BITNAMI_HOME: "${BITNAMI_HOME}"
		Description: Root directory where the Makefile resides
	
	NEO4J_HOME: "${NEO4J_HOME}"
		Description: Directory where Neo4j resides

	Common usage:

	...::: Run targets :::...
	$ make <target> <args>

endef