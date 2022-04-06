SERVICE := pipeline-functions

target:
	$(info ${HELP_MESSAGE})
	@exit 0

deploy: 
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deploying ${SERVICE}" 2>&1 | tee -a $$CFN_LOG_PATH
	$(MAKE) deploy.trigger

deploy.trigger:
	$(MAKE) -C trigger/ deploy

delete:
	$(MAKE) delete.trigger

delete.trigger:
	@echo "$$(gdate -u +'%Y-%m-%d %H:%M:%S.%3N') - Deleting ${SERVICE}" 2>&1 | tee -a $$CFN_LOG_PATH
	$(MAKE) -C trigger/ delete

#############
#  Helpers  #
#############

define HELP_MESSAGE

	Common usage:

	...::: Deploy all CloudFormation based services :::...
	$ make deploy

	...::: Delete all CloudFormation based services :::...
	$ make delete

endef