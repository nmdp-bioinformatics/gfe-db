# This Makefile only orchestrates process that are run on the EC2 database instance, it is deployed and called on the server
# Must be run as root
# TODO make sure commands are idempotent
SERVICE := database

# Application specific environment variables
# env.sh is created by the EC2 User Data in CloudFormation template
include env.sh
export

export LOGS_DIR := /tmp/logs
export REPORTS_DIR := ${EC2_USER_HOME}/reports
export NEO4J_USERNAME=$(shell aws secretsmanager get-secret-value --secret-id "/${APP_NAME}/${STAGE}/${AWS_REGION}/Neo4jCredentials" | jq -r '.SecretString' | jq -r '.NEO4J_USERNAME')
export NEO4J_PASSWORD=$(shell aws secretsmanager get-secret-value --secret-id "/${APP_NAME}/${STAGE}/${AWS_REGION}/Neo4jCredentials" | jq -r '.SecretString' | jq -r '.NEO4J_PASSWORD')

REQUIRED_VARS := APP_NAME STAGE AWS_REGION SERVICE EC2_USER_HOME NEO4J_HOME INSTANCE_ID DATA_BUCKET_NAME NEO4J_USERNAME NEO4J_PASSWORD NEO4J_DATABASE_NAME ADMIN_EMAIL GDS_VERSION NEO4J_URI NEO4J_ENCRYPTION

target:
	$(info ${HELP_MESSAGE})
	@exit 0

env:
	@printenv

get.data: #=> Download the build data locally; get.data release=3470
	# should be [[ -n "$$release" ]]
	@[ -z "$$release" ] && aws s3 cp --recursive s3://${DATA_BUCKET_NAME}/data/$$release/csv/ ${NEO4J_HOME}/import/ || echo "No release argument"

env.validate.vars:
	$(foreach var,$(REQUIRED_VARS),\
		$(if $(value $(var)),,$(error $(var) is not set. Please add $(var) to the environment variables.)))
	@echo "$$(date -u +'%Y-%m-%d %H:%M:%S.%3N') - Found environment variables" 2>&1

bootstrap: env.validate.vars
ifeq ($(USE_PRIVATE_SUBNET),false)
	$(MAKE) eip.waiter
	$(MAKE) ssl.create-cert
endif

neo4j: bootstrap
	$(MAKE) neo4j.plugins.install-gds
	$(MAKE) neo4j.restart
	$(MAKE) neo4j.init
	$(MAKE) copy-logs

# eip
eip.waiter:
	@bash ${EC2_USER_HOME}/init/eip_assoc_waiter.sh

neo4j.credentials.show:
	@echo "Username: ${NEO4J_USERNAME}"
	@echo "Password: ${NEO4J_PASSWORD}"

neo4j.status:
	@systemctl status neo4j

neo4j.stop:
	@systemctl stop neo4j

neo4j.start:
	@systemctl start neo4j
	$(MAKE) neo4j.waiter

neo4j.restart:
	@systemctl restart neo4j
	$(MAKE) neo4j.waiter

neo4j.config.backup:
	@echo "Backing up Neo4j configuration to ${NEO4J_HOME}/conf/$$(date -u +'%Y-%m-%d-%H-%M')-neo4j.conf.bkp"
	@cp "${NEO4J_HOME}/conf/neo4j.conf" "${NEO4J_HOME}/conf/$$(date -u +'%Y-%m-%d-%H-%M')-neo4j.conf.bkp"

# https://github.com/neo4j-contrib/neo4j-apoc-procedures/releases/download/5.13.0/apoc-5.13.0-extended.jar
neo4j.plugins.install-apoc-extended:
	@rm -f ${NEO4J_HOME}/plugins/apoc-*.jar
	@echo "Downloading APOC Extended plugin..."
	@curl -L https://github.com/neo4j-contrib/neo4j-apoc-procedures/releases/download/${APOC_VERSION}/apoc-${APOC_VERSION}-extended.jar -O
	@mv apoc-${APOC_VERSION}-extended.jar ${NEO4J_HOME}/plugins/apoc-${APOC_VERSION}-extended.jar

# https://graphdatascience.ninja/neo4j-graph-data-science-2.5.5.zip
neo4j.plugins.install-gds:
	@rm -f ${NEO4J_HOME}/plugins/neo4j-graph-data-science-*.jar
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
	@sudo -u neo4j mkdir -p ${NEO4J_HOME}/backups ${NEO4J_HOME}/cypher
	@echo "Fetching Cypher scripts from S3..."
	@aws s3 cp --recursive s3://${DATA_BUCKET_NAME}/${CONFIG_S3_PATH}/${SERVICE}/neo4j/cypher/ ${NEO4J_HOME}/cypher/

neo4j.query.init: # public: neo4j+s://${SUBDOMAIN}.${HOST_DOMAIN}:7687; private: bolt://<private_ip>:7687
	@echo "Executing initialization queries"
ifeq ($(USE_PRIVATE_SUBNET),true)
	$(MAKE) neo4j.query.init.no-ssl
else ifeq ($(USE_PRIVATE_SUBNET),false)
	$(MAKE) neo4j.query.init.ssl
endif

neo4j.query.init.ssl:
	@cat ${NEO4J_HOME}/cypher/init.cyp | ${NEO4J_HOME}/bin/cypher-shell -u ${NEO4J_USERNAME} -p ${NEO4J_PASSWORD} -a neo4j+s://${SUBDOMAIN}.${HOST_DOMAIN}:7687

neo4j.query.init.no-ssl:
	@cat ${NEO4J_HOME}/cypher/init.cyp | ${NEO4J_HOME}/bin/cypher-shell -u ${NEO4J_USERNAME} -p ${NEO4J_PASSWORD} -a bolt://127.0.0.1:7687 --encryption false

neo4j.query.create-users: #=> users=<username:password,...>
	@echo "Creating users"
	@bash init/create_users.sh "$$users"

neo4j.query: #=> query=<cypher string>
ifeq ($(USE_PRIVATE_SUBNET),true)
	@echo "$$query" | ${NEO4J_HOME}/bin/cypher-shell -u ${NEO4J_USERNAME} -p ${NEO4J_PASSWORD} -a bolt://127.0.0.1:7687 --encryption false
else ifeq ($(USE_PRIVATE_SUBNET),false)
	@echo "$$query" | ${NEO4J_HOME}/bin/cypher-shell -u ${NEO4J_USERNAME} -p ${NEO4J_PASSWORD} -a neo4j+s://${SUBDOMAIN}.${HOST_DOMAIN}:7687
endif

neo4j.query.show-constraints: # public: neo4j+s://${SUBDOMAIN}.${HOST_DOMAIN}:7687; private: bolt://<private_ip>:7687
ifeq ($(USE_PRIVATE_SUBNET),true)
	@echo "SHOW CONSTRAINTS;" | ${NEO4J_HOME}/bin/cypher-shell -u ${NEO4J_USERNAME} -p ${NEO4J_PASSWORD} -a bolt://127.0.0.1:7687 --encryption false
else ifeq ($(USE_PRIVATE_SUBNET),false)
	@echo "SHOW CONSTRAINTS;" | ${NEO4J_HOME}/bin/cypher-shell -u ${NEO4J_USERNAME} -p ${NEO4J_PASSWORD} -a neo4j+s://${SUBDOMAIN}.${HOST_DOMAIN}:7687
endif

neo4j.waiter:
ifeq ($(USE_PRIVATE_SUBNET),true)
	$(MAKE) neo4j.waiter.no-ssl
else ifeq ($(USE_PRIVATE_SUBNET),false)
	$(MAKE) neo4j.waiter.ssl
endif

neo4j.waiter.ssl:
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
	@cd ${NEO4J_HOME} && \
	sudo -u neo4j mkdir -p ${NEO4J_HOME}/backups/system ${NEO4J_HOME}/backups/${NEO4J_DATABASE_NAME} && \
	sudo -u neo4j ${NEO4J_HOME}/bin/neo4j-admin database dump system --to-path=${NEO4J_HOME}/backups/system/ --overwrite-destination=true --verbose && \
	sudo -u neo4j ${NEO4J_HOME}/bin/neo4j-admin database dump ${NEO4J_DATABASE_NAME} --to-path=${NEO4J_HOME}/backups/${NEO4J_DATABASE_NAME}/ --overwrite-destination=true --verbose && \
	archive_path=$$(date +%Y/%m/%d/%H/%M)/${NEO4J_DATABASE_NAME}.zip && \
	echo "Uploading archive to S3: s3://${DATA_BUCKET_NAME}/backups/${SERVICE}/${NEO4J_DATABASE_NAME}/$$archive_path" && \
	sudo -u neo4j zip -r - backups | aws s3 cp - s3://${DATA_BUCKET_NAME}/backups/${SERVICE}/${NEO4J_DATABASE_NAME}/$$archive_path && \
	echo "Updating SSM Parameter for current backup S3 path" && \
	echo "AWS region: ${AWS_REGION}" && \
	res=$$(aws ssm put-parameter \
		--name "/${APP_NAME}/${STAGE}/${AWS_REGION}/CurrentBackupS3Path" \
		--value "s3://${DATA_BUCKET_NAME}/backups/${SERVICE}/${NEO4J_DATABASE_NAME}/$$archive_path" \
		--type String \
		--overwrite \
		--region ${AWS_REGION}) && \
	echo $$res && \
	echo "Cleaning up" && \
	sudo -u neo4j rm -rf ${NEO4J_HOME}/backups/*
	$(MAKE) neo4j.start

neo4j.restore: #from_path=s3://<backup path>.zip
	$(MAKE) neo4j.stop
	@echo "Received data restore request for backup $$from_path"
	@[ "$$from_path" != "" ] || (echo "from_path is required" && exit 1)
	@echo "$$from_path" | grep -E '^s3:\/\/[a-zA-Z0-9\.-]+\/([a-zA-Z0-9\.-]+\/)*[a-zA-Z0-9\.-]+\.zip$$' > /dev/null || (echo "\`from_path\` must be an S3 URI: 's3://<bucket>/path'" && exit 1)
	@echo "Running restore for $$from_path"
	@mkdir -p /tmp/restore
	@aws s3 cp $$from_path /tmp/restore/${NEO4J_DATABASE_NAME}.zip && chown neo4j:neo4j /tmp/restore/${NEO4J_DATABASE_NAME}.zip
	@sudo -u neo4j unzip -o /tmp/restore/${NEO4J_DATABASE_NAME}.zip -d ${NEO4J_HOME}/
	@sudo -u neo4j ${NEO4J_HOME}/bin/neo4j-admin database load --from-path=${NEO4J_HOME}/backups/system system --overwrite-destination=true
	@sudo -u neo4j ${NEO4J_HOME}/bin/neo4j-admin database load --from-path=${NEO4J_HOME}/backups/${NEO4J_DATABASE_NAME} ${NEO4J_DATABASE_NAME} --overwrite-destination=true
	@rm -r /tmp/restore
	$(MAKE) neo4j.start

# # TODO update for Neo4j 5
# neo4j.create-report:
# 	@mkdir -p ${REPORTS_DIR}
# 	@${NEO4J_HOME}/bin/neo4j-admin check-consistency --database=neo4j --verbose=true --report-dir=reports > reports/consistency_check.txt

ssl.create-cert:
	@bash init/create_cert.sh "${SUBDOMAIN}.${HOST_DOMAIN}" ${ADMIN_EMAIL}

ssl.renew-cert:
	@bash init/renew_cert.sh "${SUBDOMAIN}.${HOST_DOMAIN}"

copy-logs:
	@mkdir -p ${LOGS_DIR}/var/logs
	$(MAKE) logs.bootstrap
	$(MAKE) logs.system
	@echo "Copying logs to S3..." && \
	aws s3 cp --recursive /tmp/logs/ s3://${DATA_BUCKET_NAME}/logs/database/${INSTANCE_ID}/$$(date +'%Y/%m/%d/%H')/

logs.bootstrap:
	@journalctl -b > ${LOGS_DIR}/journalctl.log

logs.system:
	@cp -r /var/log ${LOGS_DIR}/var/log

logs.bootstrap.show:
	@journalctl -b | less

logs.neo4j.show:
	@cat ${NEO4J_HOME}/logs/neo4j.log

logs.neo4j.debug.show:
	@cat ${NEO4J_HOME}/logs/debug.log

define HELP_MESSAGE

	Environment variables:

	EC2_USER_HOME: "${EC2_USER_HOME}"
		Description: Root directory where the Makefile resides
	
	NEO4J_HOME: "${NEO4J_HOME}"
		Description: Directory where Neo4j resides

	Common usage:

	...::: Run targets :::...
	$ make <target> <args>

endef