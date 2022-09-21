# This Makefile only orchestrates process that are run on the EC2 database instance, it is deployed and called on the server
# Run as sudo only

# Application specific environment variables
include env.sh
export

export ROOT_DIR=/home/bitnami
export NEO4J_HOME=/opt/bitnami/neo4j
export NEO4J_USERNAME=$(shell cat ${ROOT_DIR}/bitnami_credentials | sed 's/ /\n/g' | grep -E "'([a-zA-Z0-9]+)'" | sed "s/'//g" | sed "s/\.//g" | head -n 1)
export NEO4J_PASSWORD=$(shell cat ${ROOT_DIR}/bitnami_credentials | sed 's/ /\n/g' | grep -E "'([a-zA-Z0-9]+)'" | sed "s/'//g" | sed "s/\.//g" | tail -n 1)

target:
	$(info ${HELP_MESSAGE})
	@exit 0

env:
	@printenv

# TODO fetch all shell scripts from S3
# get.scripts:
# get.cypher:
get.data: #=> Download the build data locally; get.data release=3470
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
	$(MAKE) neo4j.publish-credentials
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
	@bash ${ROOT_DIR}/init/eip_assoc_waiter.sh

# secretsmanager
neo4j.publish-credentials:
	@aws --region ${AWS_REGION} secretsmanager update-secret \
		--secret-id ${NEO4J_CREDENTIALS_SECRET_ARN} \
		--secret-string "{\"NEO4J_USERNAME\":\"${NEO4J_USERNAME}\",\"NEO4J_PASSWORD\":\"${NEO4J_PASSWORD}\"}"

# neo4j
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
	@timeout=300 && \
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
	@cd ${BITNAMI_NEO4J}/data/ && zip ${BITNAMI_HOME}/gfedb.zip databases/
	@result=$$(aws s3 cp ${BITNAMI_HOME}/gfedb.zip s3://${DATA_BUCKET_NAME}/backups/neo4j/$$(date +'%Y/%m/%d/%H')/gfedb.zip) && \
		echo $$result
	@rm ${BITNAMI_HOME}/gfedb.zip

# neo4j.restore: #s3_path
# # get path (path selection process should get most recent)
# # download archive from s3
# # unzip and overwrite
# 	@unzip -o ${BITNAMI_HOME}/gfedb.zip 
# # -d data/databases/

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
	@echo "Copying logs to S3..."
	@mkdir -p /tmp/logs
	@export instance_id=$$(curl http://169.254.169.254/latest/meta-data/instance-id) && \
	journalctl -b > /tmp/logs/$$instance_id-boot.log
	@echo "Copying Neo4j logs to S3..."
	@journalctl -e -u neo4j > /tmp/logs/neo4j.log
	@aws s3 cp --recursive /tmp/logs/ s3://${DATA_BUCKET_NAME}/logs/database/bootstrap/$$(date +'%Y/%m/%d/%H')/

define HELP_MESSAGE

	Environment variables:

	ROOT_DIR: "${ROOT_DIR}"
		Description: Root directory where the Makefile resides
	
	NEO4J_HOME: "${NEO4J_HOME}"
		Description: Directory where Neo4j resides

	Common usage:

	...::: Run targets :::...
	$ make <target> <args>

endef