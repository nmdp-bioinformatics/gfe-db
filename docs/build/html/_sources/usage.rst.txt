Usage
=====

Deployment to AWS
-----------------

Follow the steps to build and deploy the application to AWS.

Quick Start
~~~~~~~~~~~

This list outlines the basic steps for deployment. For more details
please see the following sections.

#. Install `prerequisites <prerequisites_>`__ 
#. Set `environment variables <environment_>`__ 
#. Check the config JSONs (parameters and state) and edit the values as desired 
#. Run ``make deploy`` to deploy the stacks to AWS 
#. Run ``make load.database release=<version>`` to load the Neo4j 
#. Run ``make get.neo4j`` to get the URL for the Neo4j browser

.. _prerequisites:

Prerequisites
~~~~~~~~~~~~~

Please refer to the respective documentation for specific installation
instructions.

* GNU Make 3.81
* coreutils (optional)
* AWS CLI
* SAM CLI
* Docker
* jq

.. _environment:

Environment
~~~~~~~~~~~

AWS Credentials
^^^^^^^^^^^^^^^

Valid AWS credentials must be available to AWS CLI and SAM CLI. The
easiest way to do this is with the following steps. 
#. Run ``aws configure`` and follow the prompts, or copy/paste them into ``~/.aws/credentials`` 
#. Export the ``AWS_PROFILE`` variable for the chosen profile to the shell environment.

.. code:: bash

   export AWS_PROFILE=default

For more information visit the documentation page: `Configuration and credential file settings 
<https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-files.html>`__

Shell Variables
^^^^^^^^^^^^^^^

These variables must be present in the shell environment before running
Make. The best way to set these variables is with a ``.env`` file
following these steps.

1. Create a ``.env`` file in the project root and add the values.

.. code:: bash

   STAGE=<dev or prod>
   APP_NAME=gfe-db
   REGION=<AWS region>
   NEO4J_USERNAME=<secret>
   NEO4J_PASSWORD=<secret>
   GITHUB_PERSONAL_ACCESS_TOKEN=<secret>

2. Source the variables to the environment.

.. code:: bash

   # Export .env file variables
   set -a && source .env && set +a

   # Check that the variables were set
   env

.. important::

    Always use a ``.env`` file, AWS SSM Parameter Store, or
    Secrets Manager for sensitive variables like credentials and API keys.
    Never hard-code them, including when developing. AWS will quarantine an
    account if any credentials get accidentally exposed and this will cause
    problems. Make sure to update ``.gitignore`` to avoid pushing sensitive
    data to public repositories.

Makefile Usage
~~~~~~~~~~~~~~

Once an AWS profile is configured and environment variables are
exported, the application can be deployed using ``make``.

.. code:: bash

   make deploy

It is also possible to deploy or update individual services.

.. code:: bash

   # Deploy/update only the infrastructure service
   make deploy.infrastructure

   # Deploy/update only the database service
   make deploy.database

   # Deploy/update only the pipeline service
   make deploy.pipeline

.. note::
    It is recommended to only deploy from the project root. This is
    because common parameters are passed from the root Makefile to nested
    Makefiles. If a stack has not been changed, the deployment script will
    continue until it reaches a stack with changes and deploy that.

Makefile Command Reference
^^^^^^^^^^^^^^^^^^^^^^^^^^

To see a list of possible commands using Make, run ``make`` on the
command line.

.. code:: bash

   # Deploy all CloudFormation based services
   make deploy

   # Deploy config files and scripts to S3
   make deploy.config

   # Run the StepFunctions State Machine to load Neo4j
   make load.database releases=<version> align=<boolean> kir=<boolean> limit=<int>

   # Download CSV data from S3
   make get.data

   # Download logs from EC2
   make get.logs

   # Display the Neo4j Browser endpoint URL
   make get.neo4j

   # Delete all CloudFormation based services and data
   make delete
