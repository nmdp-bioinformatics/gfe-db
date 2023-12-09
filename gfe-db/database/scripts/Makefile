# This Makefile only orchestrates process that are run on the EC2 database instance, it is deployed and called on the server
# Must be run as sudo

# Application specific environment variables
include env.sh
export

# export BITNAMI_HOME ?= /home/bitnami
export LOGS_DIR := /tmp/logs
export REPORTS_DIR := ${BITNAMI_HOME}/reports
export NEO4J_HOME := /opt/bitnami/neo4j
export NEO4J_USERNAME=$(shell awk -F"'" '/The default username and password is/ {print $$2}' bitnami_credentials)
export NEO4J_PASSWORD=$(shell awk -F"'" '/The default username and password is/ {print $$4}' bitnami_credentials)
# export INSTANCE_ID := $(shell curl -sS http://169.254.169.254/latest/meta-data/instance-id)

REQUIRED_VARS := BITNAMI_HOME NEO4J_HOME INSTANCE_ID DATA_BUCKET_NAME NEO4J_USERNAME NEO4J_PASSWORD ADMIN_EMAIL GDS_VERSION

target:
	$(info ${HELP_MESSAGE})
	@exit 0

env:
	@printenv

get.data: #=> Download the build data locally; get.data release=3470
	# should be [[ -n "$$release" ]]
	@[ -z "$$release" ] && aws s3 cp --recursive s3://${DATA_BUCKET_NAME}/data/$$release/csv/ ${NEO4J_HOME}/import/ || echo "No release argument"
# env.fetch:

env.validate.vars:
	$(foreach var,$(REQUIRED_VARS),\
		$(if $(value $(var)),,$(error $(var) is not set. Please add $(var) to the environment variables.)))
	@echo "$$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Found environment variables" 2>&1

# env.check:
# ifndef DATA_BUCKET_NAME
# $(error DATA_BUCKET_NAME is not set.)
# endif
# ifndef NEO4J_HOME
# $(error NEO4J_HOME is not set.)
# endif
# ifndef NEO4J_USERNAME
# $(error NEO4J_HOME is not set.)
# endif
# ifndef NEO4J_PASSWORD
# $(error NEO4J_HOME is not set.)
# endif
# # ifndef HOST_DOMAIN
# # $(error HOST_DOMAIN is not set.)
# # endif
# ifndef ADMIN_EMAIL
# $(error ADMIN_EMAIL is not set.)
# endif
# ifndef APOC_VERSION
# $(error APOC_VERSION is not set.)
# endif
# ifndef GDS_VERSION
# $(error GDS_VERSION is not set.)
# endif
# 	@echo "$$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Found environment variables" 2>&1

bootstrap: env.validate.vars
	# $(MAKE) cfn-helpers.install
	# $(MAKE) ssm.install
	# $(MAKE) cloudwatch.install
	# $(MAKE) cfn-signal exitcode=0
ifeq ($(USE_PRIVATE_SUBNET),false)
	$(MAKE) eip.waiter
	$(MAKE) ssl.create-cert
endif

neo4j: bootstrap
	$(MAKE) neo4j.update-credentials
	$(MAKE) neo4j.config.update
	$(MAKE) neo4j.plugins.install-gds
	$(MAKE) neo4j.restart
	$(MAKE) neo4j.init
	$(MAKE) neo4j.restart
	$(MAKE) copy-logs

# eip
eip.waiter:
	@bash ${BITNAMI_HOME}/init/eip_assoc_waiter.sh

# neo4j
neo4j.update-credentials:
	@aws --region ${AWS_REGION} secretsmanager update-secret \
		--secret-id "/$${APP_NAME}/$${STAGE}/$${AWS_REGION}/Neo4jCredentials" \
		--secret-string "{\"NEO4J_USERNAME\":\"${NEO4J_USERNAME}\",\"NEO4J_PASSWORD\":\"${NEO4J_PASSWORD}\"}"

neo4j.credentials.show:
	@echo "Username: ${NEO4J_USERNAME}"
	@echo "Password: ${NEO4J_PASSWORD}"

neo4j.status:
	@/opt/bitnami/ctlscript.sh status

neo4j.stop:
	@/opt/bitnami/ctlscript.sh stop

neo4j.start:
	$(MAKE) neo4j.update-permissions
	@/opt/bitnami/ctlscript.sh start
	$(MAKE) neo4j.waiter

neo4j.restart:
	$(MAKE) neo4j.update-permissions
	@/opt/bitnami/ctlscript.sh restart
	$(MAKE) neo4j.waiter

neo4j.config.update:
	$(MAKE) neo4j.config.download
	$(MAKE) neo4j.config.update-dns

# TODO test init scripts for public and private deployment
neo4j.config.download: neo4j.config.backup
	@echo "Downloading Neo4j configuration from ${DATA_BUCKET_NAME}"
	@aws s3 cp s3://${DATA_BUCKET_NAME}/config/neo4j/neo4j.conf ${NEO4J_HOME}/conf/neo4j.template

# TODO test init scripts for public and private deployment
neo4j.config.update-dns:
ifeq ($(USE_PRIVATE_SUBNET),true)
	$(info Updating Neo4j config to use private DNS)
	$(MAKE) neo4j.config.update-dns.private
else ifeq ($(USE_PRIVATE_SUBNET),false)
	$(info Updating Neo4j config to use public DNS)
	$(MAKE) neo4j.config.update-dns.public
endif

# This target will configure the database to use the correct DNS name using env vars
neo4j.config.update-dns.private:
	@cat ${NEO4J_HOME}/conf/neo4j.template | \
	sed "s|# server.default_advertised_address=127\.0\.0\.1|server.default_advertised_address=${PRIVATE_IP}|" > ${NEO4J_HOME}/conf/neo4j.conf

# This target will configure the database to use the correct DNS name using env vars
neo4j.config.update-dns.public:
	@[ "${HOST_DOMAIN}" != "" ] && \
	cat ${NEO4J_HOME}/conf/neo4j.template | \
	sed s"|# server.default_advertised_address=127\.0\.0\.1|server.default_advertised_address=${SUBDOMAIN}.${HOST_DOMAIN}|" > ${NEO4J_HOME}/conf/neo4j.conf

# TODO update for Neo4j 5
neo4j.config.reset-auth.ssl:
	@echo "Resetting Neo4j authentication"
	$(MAKE) neo4j.stop
	@sed -i 's/#dbms.security.auth_enabled=false/dbms.security.auth_enabled=false/g' ${NEO4J_HOME}/conf/neo4j.conf
	$(MAKE) neo4j.start
	@-echo "ALTER USER neo4j SET PASSWORD '${NEO4J_PASSWORD}'" | /opt/bitnami/neo4j/bin/cypher-shell -u neo4j -p neo4j -a neo4j+s://${SUBDOMAIN}.${HOST_DOMAIN}:7687 -d system
	$(MAKE) neo4j.stop
	@sed -i 's/dbms.security.auth_enabled=false/#dbms.security.auth_enabled=false/g' ${NEO4J_HOME}/conf/neo4j.conf
	$(MAKE) neo4j.start

# TODO update for Neo4j 5
neo4j.config.reset-auth.no-ssl:
	@echo "Resetting Neo4j authentication"
	$(MAKE) neo4j.stop
	@sed -i 's/#dbms.security.auth_enabled=false/dbms.security.auth_enabled=false/g' ${NEO4J_HOME}/conf/neo4j.conf
	$(MAKE) neo4j.start
	@-echo "ALTER USER neo4j SET PASSWORD '${NEO4J_PASSWORD}'" | /opt/bitnami/neo4j/bin/cypher-shell -u neo4j -p neo4j -a bolt://127.0.0.1:7687 --encryption false -d system
	$(MAKE) neo4j.stop
	@sed -i 's/dbms.security.auth_enabled=false/#dbms.security.auth_enabled=false/g' ${NEO4J_HOME}/conf/neo4j.conf
	$(MAKE) neo4j.start

neo4j.config.backup:
	@echo "Backing up Neo4j configuration to ${NEO4J_HOME}/conf/$$(date -u +'%Y-%m-%d-%H-%M')-neo4j.conf.bkp"
	@cp "${NEO4J_HOME}/conf/neo4j.conf" "${NEO4J_HOME}/conf/$$(date -u +'%Y-%m-%d-%H-%M')-neo4j.conf.bkp"

# TODO check if the APOC plugin is already installed for the version specified and skip if it is
# https://github.com/neo4j-contrib/neo4j-apoc-procedures/releases/download/5.13.0/apoc-5.13.0-extended.jar
neo4j.plugins.install-apoc-extended:
	@rm -f ${NEO4J_HOME}/plugins/apoc-*.jar
	@echo "Downloading APOC Extended plugin..."
	@curl -L https://github.com/neo4j-contrib/neo4j-apoc-procedures/releases/download/${APOC_VERSION}/apoc-${APOC_VERSION}-extended.jar -O
	@mv apoc-${APOC_VERSION}-extended.jar ${NEO4J_HOME}/plugins/apoc-${APOC_VERSION}-extended.jar

# TODO check if the GDS plugin is already installed for the version specified and skip if it is
# https://graphdatascience.ninja/neo4j-graph-data-science-2.5.5.zip
neo4j.plugins.install-gds:
	@echo "Downloading Neo4j Graph Data Science plugin..."
	@curl -L https://graphdatascience.ninja/neo4j-graph-data-science-${GDS_VERSION}.zip -O
	@unzip neo4j-graph-data-science-${GDS_VERSION}.zip
	@mv neo4j-graph-data-science-${GDS_VERSION}.jar ${NEO4J_HOME}/plugins/neo4j-graph-data-science-${GDS_VERSION}.jar
	@echo "Cleaning up"
	@rm neo4j-graph-data-science-${GDS_VERSION}.zip

neo4j.init:
	$(MAKE) neo4j.init.download-scripts
	$(MAKE) neo4j.query.init

neo4j.update-permissions:
	@chown -R neo4j:neo4j ${NEO4J_HOME}/*

neo4j.init.download-scripts:
	@mkdir -p ${NEO4J_HOME}/backups ${NEO4J_HOME}/cypher
	@echo "Fetching Cypher scripts from S3..."
	@aws s3 cp --recursive s3://${DATA_BUCKET_NAME}/config/neo4j/cypher/ ${NEO4J_HOME}/cypher/

neo4j.query.init: # public: neo4j+s://${SUBDOMAIN}.${HOST_DOMAIN}:7687; private: bolt://<private_ip>:7687
ifeq ($(USE_PRIVATE_SUBNET),true)
	$(MAKE) neo4j.query.init.no-ssl
else ifeq ($(USE_PRIVATE_SUBNET),false)
	$(MAKE) neo4j.query.init.ssl
endif

neo4j.query.init.ssl:
	@echo "Executing initialization queries..."
	@cat ${NEO4J_HOME}/cypher/init.cyp | /opt/bitnami/neo4j/bin/cypher-shell -u ${NEO4J_USERNAME} -p ${NEO4J_PASSWORD} -a neo4j+s://${SUBDOMAIN}.${HOST_DOMAIN}:7687

neo4j.query.show-constraints: # public: neo4j+s://${SUBDOMAIN}.${HOST_DOMAIN}:7687; private: bolt://<private_ip>:7687
ifeq ($(USE_PRIVATE_SUBNET),true)
	@echo "SHOW CONSTRAINTS;" | /opt/bitnami/neo4j/bin/cypher-shell -u ${NEO4J_USERNAME} -p ${NEO4J_PASSWORD} -a bolt://127.0.0.1:7687 --encryption false
else ifeq ($(USE_PRIVATE_SUBNET),false)
	@echo "SHOW CONSTRAINTS;" | /opt/bitnami/neo4j/bin/cypher-shell -u ${NEO4J_USERNAME} -p ${NEO4J_PASSWORD} -a neo4j+s://${SUBDOMAIN}.${HOST_DOMAIN}:7687
endif

neo4j.query.init.no-ssl:
	@echo "Executing initialization queries..."
	@cat ${NEO4J_HOME}/cypher/init.cyp | /opt/bitnami/neo4j/bin/cypher-shell -u ${NEO4J_USERNAME} -p ${NEO4J_PASSWORD} -a bolt://127.0.0.1:7687 --encryption false

neo4j.waiter:
ifeq ($(USE_PRIVATE_SUBNET),true)
	$(MAKE) neo4j.waiter.no-ssl
else ifeq ($(USE_PRIVATE_SUBNET),false)
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
endif

neo4j.waiter.no-ssl:
	@timeout=120 && \
	counter=0 && \
	host=127.0.0.1 && \
	echo "Waiting for response from Neo4j at http://$$host:7474..." && \
	until $$(curl --output /dev/null --silent --head --fail http://$$host:7474) ; do \
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

# TODO update for Neo4j 5
# TODO testing s3://prod-gfe-db-531868584498-us-east-1/backups/neo4j/2023/09/14/05
neo4j.restore: #from_path=s3://<backup path>
	@echo "Received data restore request for backup $$from_path"
	@[ "$$from_path" != "" ] || (echo "from_path is required" && exit 1)
	@echo "$$from_path" | grep -E '^s3:\/\/[a-zA-Z0-9\.-]+\/([a-zA-Z0-9\.-]+\/)*[a-zA-Z0-9\.-]+\.zip$$' > /dev/null || (echo "\`from_path\` must be an S3 URI: 's3://<bucket>/path'" && exit 1)
	@echo "Running restore for $$from_path"
	$(MAKE) neo4j.stop
	@aws s3 cp $$from_path ${BITNAMI_HOME}/gfedb.zip
	@unzip -o ${BITNAMI_HOME}/gfedb.zip
	@cp -r data/* ${BITNAMI_NEO4J}/data/
	@chown -R neo4j:neo4j ${BITNAMI_NEO4J}/data/
	@rm ${BITNAMI_HOME}/gfedb.zip
	@rm -r ${BITNAMI_HOME}/data/
ifeq ($(USE_PRIVATE_SUBNET),true)
	$(MAKE) neo4j.config.reset-auth.no-ssl
else ifeq ($(USE_PRIVATE_SUBNET),false)
	$(MAKE) neo4j.config.reset-auth.ssl
endif

neo4j.query.create-constraints:
	@mkdir -p ${NEO4J_HOME}/backups ${NEO4J_HOME}/cypher
	@echo "Fetching cypher scripts from S3..."
	@aws s3 cp --recursive s3://${DATA_BUCKET_NAME}/config/neo4j/cypher/ ${NEO4J_HOME}/cypher/
	@cat ${NEO4J_HOME}/cypher/create_constraints.cyp | cypher-shell -u ${NEO4J_USERNAME} -p ${NEO4J_PASSWORD} -a neo4j+s://${SUBDOMAIN}.${HOST_DOMAIN}:7687

neo4j.query.drop-constraints:
	@mkdir -p ${NEO4J_HOME}/backups ${NEO4J_HOME}/cypher
	@echo "Fetching cypher scripts from S3..."
	@aws s3 cp --recursive s3://${DATA_BUCKET_NAME}/config/neo4j/cypher/ ${NEO4J_HOME}/cypher/
	@cat ${NEO4J_HOME}/cypher/drop_constraints.cyp | cypher-shell -u ${NEO4J_USERNAME} -p ${NEO4J_PASSWORD} -a neo4j+s://${SUBDOMAIN}.${HOST_DOMAIN}:7687

# TODO update for Neo4j 5
neo4j.create-report:
	@mkdir -p ${REPORTS_DIR}
	@neo4j-admin check-consistency --database=neo4j --verbose=true --report-dir=reports > reports/consistency_check.txt

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
	@cp -r /var/log/* ${LOGS_DIR}/system/

logs.neo4j:
	@cp -r ${NEO4J_HOME}/logs/* ${LOGS_DIR}/neo4j/

logs.bootstrap.show:
	@journalctl -b | less

logs.neo4j.show:
	@cat ${NEO4J_HOME}/logs/neo4j.log

logs.neo4j.debug.show:
	@cat ${NEO4J_HOME}/logs/debug.log

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