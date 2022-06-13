# This Makefile only orchestrates process that are run on the EC2 database instance, it is deployed and called on the server
# Run as sudo only

# Application specific environment variables
include .env
export

export ROOT_DIR=/home/bitnami
export NEO4J_HOME=/opt/bitnami/neo4j
export NEO4J_USERNAME=$(shell cat ${ROOT_DIR}/bitnami_credentials | sed 's/ /\n/g' | grep -E "'([a-zA-Z0-9]+)'" | sed "s/'//g" | sed "s/\.//g" | head -n 1)
export NEO4J_PASSWORD=$(shell cat ${ROOT_DIR}/bitnami_credentials | sed 's/ /\n/g' | grep -E "'([a-zA-Z0-9]+)'" | sed "s/'//g" | sed "s/\.//g" | tail -n 1)

target:
	$(info ${HELP_MESSAGE})
	@exit 0

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
	$(MAKE) eip.waiter
	$(MAKE) ssl.create-cert

neo4j: bootstrap
	$(MAKE) neo4j.publish-credentials
	$(MAKE) neo4j.config.update
	$(MAKE) neo4j.plugins.install-apoc
	$(MAKE) neo4j.plugins.install-gds
	$(MAKE) neo4j.restart
	$(MAKE) neo4j.init
	$(MAKE) copy-logs

# eip
eip.waiter:
	@bash ${ROOT_DIR}/init/eip_assoc_waiter.sh

# secretsmanager
neo4j.publish-credentials:
	@aws --region ${REGION} secretsmanager update-secret \
		--secret-id ${NEO4J_CREDENTIALS_SECRET_ARN} \
		--secret-string "{\"NEO4J_USERNAME\":\"${NEO4J_USERNAME}\",\"NEO4J_PASSWORD\":\"${NEO4J_PASSWORD}\"}"

# neo4j
neo4j.status:
	@/opt/bitnami/ctlscript.sh status

neo4j.stop:
	@/opt/bitnami/ctlscript.sh stop

neo4j.start:
	@/opt/bitnami/ctlscript.sh start
	@timeout=300 && \
	counter=0 && \
	until $$(curl --output /dev/null --silent --head --fail http://localhost:7474) ; do \
		printf '.' ; \
		sleep 1 ; \
		counter=$$((counter + 1)) ; \
		[ $$counter -eq $$timeout ] && exit 1 || true ; \
	done
	@printf "%s\n" " "
	@echo "Neo4j is ready"

neo4j.restart:
	@/opt/bitnami/ctlscript.sh restart
	@timeout=300 && \
	counter=0 && \
	until $$(curl --output /dev/null --silent --head --fail http://localhost:7474) ; do \
		printf '.' ; \
		sleep 1 ; \
		counter=$$((counter + 1)) ; \
		[ $$counter -eq $$timeout ] && exit 1 || true ; \
	done
	@printf "%s\n" " "
	@echo "Neo4j is ready"

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
	@mkdir -p ${NEO4J_HOME}/{backups,cypher}
	@echo "Fetching cypher scripts from S3..."
	@aws s3 cp --recursive s3://${DATA_BUCKET_NAME}/config/neo4j/cypher/ ${NEO4J_HOME}/cypher/
	@cat ${NEO4J_HOME}/cypher/init.cyp | /opt/bitnami/neo4j/bin/cypher-shell -u ${NEO4J_USERNAME} -p ${NEO4J_PASSWORD} -a neo4j+s://${HOST_DOMAIN}:7687

# TODO test
neo4j.import:
	@mkdir -p ${NEO4J_HOME}/backups ${NEO4J_HOME}/cypher
	@echo "Fetching cypher scripts from S3..."
	@aws s3 cp --recursive s3://${DATA_BUCKET_NAME}/config/neo4j/cypher/ ${NEO4J_HOME}/cypher/
	@cat ${NEO4J_HOME}/cypher/${Neo4jInitScript} | cypher-shell -u neo4j -p ${NEO4J_PASSWORD}

neo4j.create-constraints:
	@mkdir -p ${NEO4J_HOME}/{backups,cypher}
	@echo "Fetching cypher scripts from S3..."
	@aws s3 cp --recursive s3://${DATA_BUCKET_NAME}/config/neo4j/cypher/ ${NEO4J_HOME}/cypher/
	@cat ${NEO4J_HOME}/cypher/create_constraints.cyp | cypher-shell -u ${NEO4J_USERNAME} -p ${NEO4J_PASSWORD} -a neo4j+s://${HOST_DOMAIN}:7687

neo4j.drop-constraints:
	@mkdir -p ${NEO4J_HOME}/{backups,cypher}
	@echo "Fetching cypher scripts from S3..."
	@aws s3 cp --recursive s3://${DATA_BUCKET_NAME}/config/neo4j/cypher/ ${NEO4J_HOME}/cypher/
	@cat ${NEO4J_HOME}/cypher/drop_constraints.cyp | cypher-shell -u ${NEO4J_USERNAME} -p ${NEO4J_PASSWORD} -a neo4j+s://${HOST_DOMAIN}:7687

ssl.create-cert:
	@bash init/create_cert.sh ${HOST_DOMAIN} ${ADMIN_EMAIL}

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