target:
	$(info ${HELP_MESSAGE})
	@exit 0

deploy: deploy.trigger

deploy.trigger:
	$(MAKE) -C trigger/ deploy

delete: delete.trigger

delete.trigger:
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