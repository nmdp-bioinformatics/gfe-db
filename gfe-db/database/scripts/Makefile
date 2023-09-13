# This Makefile only orchestrates process that are run on the EC2 database instance, it is deployed and called on the server
# Must be run as sudo

# Application specific environment variables
include env.sh
export

export BITNAMI_HOME ?= /home/bitnami
export LOGS_DIR := /tmp/logs
export REPORTS_DIR := ${BITNAMI_HOME}/reports
export NEO4J_HOME := /opt/bitnami/neo4j
export NEO4J_USERNAME=$(shell awk -F"'" '/The default username and password is/ {print $$2}' bitnami_credentials)
export NEO4J_PASSWORD=$(shell awk -F"'" '/The default username and password is/ {print $$4}' bitnami_credentials)
export INSTANCE_ID := $(shell curl -sS http://169.254.169.254/latest/meta-data/instance-id)

target:
	$(info ${HELP_MESSAGE})
	@exit 0

env:
	@printenv

get.data: #=> Download the build data locally; get.data release=3470
	# should be [[ -n "$$release" ]]
	@[ -z "$$release" ] && aws s3 cp --recursive s3://${DATA_BUCKET_NAME}/data/$$release/csv/ ${NEO4J_HOME}/import/ || echo "No release argument"

env.check:
ifndef DATA_BUCKET_NAME
$(error DATA_BUCKET_NAME is not set.)
endif
ifndef NEO4J_HOME
$(error NEO4J_HOME is not set.)
endif
ifndef NEO4J_USERNAME
$(error NEO4J_USERNAME is not set.)
endif
ifndef NEO4J_PASSWORD
$(error NEO4J_PASSWORD is not set.)
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
	@echo "$$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Found environment variables" 2>&1

bootstrap:
	# $(MAKE) cfn-helpers.install
	# $(MAKE) ssm.install
	# $(MAKE) cloudwatch.install
	# $(MAKE) cfn-signal exitcode=0
	$(MAKE) eip.waiter
	$(MAKE) ssl.create-cert

neo4j: bootstrap
	$(MAKE) neo4j.credentials.update
	$(MAKE) neo4j.config.update
	# $(MAKE) neo4j.plugins.install-apoc # TODO already installed
	# $(MAKE) neo4j.plugins.install-gds # TODO enable
	$(MAKE) neo4j.restart
	$(MAKE) neo4j.init
	$(MAKE) neo4j.restart
	$(MAKE) copy-logs

# eip
eip.waiter:
	@bash ${BITNAMI_HOME}/init/eip_assoc_waiter.sh

# neo4j
neo4j.credentials.update:
	@aws --region ${AWS_REGION} secretsmanager update-secret \
		--secret-id ${NEO4J_CREDENTIALS_SECRET_ARN} \
		--secret-string "{\"NEO4J_USERNAME\":\"${NEO4J_USERNAME}\",\"NEO4J_PASSWORD\":\"${NEO4J_PASSWORD}\"}"

neo4j.credentials.show:
	@echo "Username: ${NEO4J_USERNAME}"
	@echo "Password: ${NEO4J_PASSWORD}"

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

neo4j.config.update: neo4j.config.backup
	@echo "Downloading Neo4j configuration from ${DATA_BUCKET_NAME}"
	@aws s3 cp s3://${DATA_BUCKET_NAME}/config/neo4j/neo4j.conf ${NEO4J_HOME}/conf/neo4j.conf

neo4j.config.reset-auth:
	@echo "Resetting Neo4j authentication"
	$(MAKE) neo4j.stop
	@sed -i 's/#dbms.security.auth_enabled=false/dbms.security.auth_enabled=false/g' ${NEO4J_HOME}/conf/neo4j.conf
	$(MAKE) neo4j.start
	@-echo "ALTER USER neo4j SET PASSWORD '${NEO4J_PASSWORD}'" | /opt/bitnami/neo4j/bin/cypher-shell -u neo4j -p neo4j -a neo4j+s://${SUBDOMAIN}.${HOST_DOMAIN}:7687 -d system
	$(MAKE) neo4j.stop
	@sed -i 's/dbms.security.auth_enabled=false/#dbms.security.auth_enabled=false/g' ${NEO4J_HOME}/conf/neo4j.conf
	$(MAKE) neo4j.start

neo4j.config.backup:
	@echo "Backing up Neo4j configuration to ${NEO4J_HOME}/conf/$$(date -u +'%Y-%m-%d-%H-%M')-neo4j.conf.bkp"
	@cp "${NEO4J_HOME}/conf/neo4j.conf" "${NEO4J_HOME}/conf/$$(date -u +'%Y-%m-%d-%H-%M')-neo4j.conf.bkp"

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

neo4j.backup:
	$(MAKE) neo4j.stop
	@cd ${BITNAMI_NEO4J}/ && zip -r ${BITNAMI_HOME}/gfedb.zip data
	@result=$$(aws s3 cp ${BITNAMI_HOME}/gfedb.zip s3://${DATA_BUCKET_NAME}/backups/neo4j/$$(date +'%Y/%m/%d/%H')/gfedb.zip) && \
		echo $$result
	@rm ${BITNAMI_HOME}/gfedb.zip
	$(MAKE) neo4j.start

neo4j.restore: #from_date=<YYYY/MM/DD/HH>
	@echo "Received data restore request for date $$from_date"
	@[ "$$from_date" != "" ] || (echo "from_date is required" && exit 1) && \
	[ "$$(echo "$$from_date" | grep -E '^[0-9]{4}/[0-9]{2}/[0-9]{2}/[0-9]{2}$$')" ] || (echo "from_date must be in the format YYYY/MM/DD/HH" && exit 1) && \
	echo "Running restore for $$from_date"
	$(MAKE) neo4j.stop
	@aws s3 cp s3://${DATA_BUCKET_NAME}/backups/neo4j/$$from_date/gfedb.zip ${BITNAMI_HOME}/gfedb.zip
	@unzip -o ${BITNAMI_HOME}/gfedb.zip
	@cp -r data/* ${BITNAMI_NEO4J}/data/
	@chown -R neo4j:neo4j ${BITNAMI_NEO4J}/data/
	@rm ${BITNAMI_HOME}/gfedb.zip
	@rm -r ${BITNAMI_HOME}/data/
	$(MAKE) neo4j.config.reset-auth

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

neo4j.report.check-consistency:
	$(MAKE) neo4j.stop
	@-echo "Removing previous reports..." && \
	rm -r reports/check-consistency/*.report
	@neo4j-admin database check neo4j --verbose --report-path reports/check-consistency/
	$(MAKE) neo4j.start

neo4j.copy-reports:
	@echo "Copying reports to S3..."
	@aws s3 cp --recursive reports/ s3://${DATA_BUCKET_NAME}/reports/database/${INSTANCE_ID}/$$(date +'%Y/%m/%d/%H')/

ssl.create-cert:
	@bash init/create_cert.sh "${SUBDOMAIN}.${HOST_DOMAIN}" ${ADMIN_EMAIL}

ssl.renew-cert:
	@bash init/renew_cert.sh "${SUBDOMAIN}.${HOST_DOMAIN}"

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