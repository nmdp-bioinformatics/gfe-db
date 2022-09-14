Usage
=====

Deployment to AWS
-----------------

Follow the steps to build and deploy the application to AWS. For information 
about the architecture deployed please see :ref:`architecture`.

Quick Start
~~~~~~~~~~~

This list outlines the basic steps for deployment. For more details
please see the following sections.

#. Install `prerequisites <prerequisites_>`__ 
#. Set `environment variables <environment_>`__ 
#. Check the config JSONs (parameters and state) and edit the values as desired 
#. Run ``make deploy`` to deploy the stacks to AWS 
#. Run ``make load.database release=<version>`` to load Neo4j
#. Run ``make get.neo4j`` to get the URL for the Neo4j browser

.. _prerequisites:

Prerequisites
~~~~~~~~~~~~~

Please refer to the respective documentation for specific installation
instructions.

* GNU Make 3.81
* coreutils (optional for macOS)
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
   AWS_REGION=<AWS region>
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
exported, the application can be deployed using ``make``. For more `make`
commands please see :ref:`makefileref`.

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

Deploying Configuration Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To deploy updates to state and/or pipeline input parameters, run the
command.

.. code:: bash

   make deploy.config

Clean Up
~~~~~~~~

To tear down resources run the command. You will need to manually delete
the data in the S3 bucket first.

.. code:: bash

   make delete

Use the following commands to tear down individual services.

.. code:: bash

   # Delete only the infrastructure service. 
   make delete.infrastructure

   # Delete only the database service
   make delete.database

   # Delete only the pipeline service
   make delete.pipeline

.. warning::
   Deleting and re-deploying a layer may cause the parameters shared by
   other layers to go out of date. To avoid this the recommendation is to deploy
   sequentially and teardown in reverse sequence. For example, tearing down and
   re-deploying the database stack may affect parameters shared with the pipeline stack, 
   causing the pipeline stack to fail. The solution would be to also be 
   tear down and re-deploy the pipeline stack.

Loading Releases
----------------

Input Parameters
~~~~~~~~~~~~~~~~

Base input parameters (excluding the ``releases`` value) are passed to
the Step Functions State Machine and determine it's behavior during
build. These are stored in a configuration file in S3 (see the :ref:`datapipelineconfig` 
configuration reference) but can be overridden. The ``releases`` value 
is appended at runtime by the trigger Lambda when it finds a new release 
in the source repository. 

+-------------+---------------+--------+--------------------------------------------------------------------+
| Variable    | Example Value | Type   | Description                                                        |
+=============+===============+========+====================================================================+
| LIMIT       | 100           | string | Number of alleles to build. Leave blank ("") to build all alleles. |
+-------------+---------------+--------+--------------------------------------------------------------------+
| ALIGN       | False         | string | Include or exclude alignments in the build                         |
+-------------+---------------+--------+--------------------------------------------------------------------+
| KIR         | False         | string | Include or exclude KIR data alignments in the build                |
+-------------+---------------+--------+--------------------------------------------------------------------+
| MEM_PROFILE | False         | string | Enable memory profiling (for catching memory leaks during build)   |
+-------------+---------------+--------+--------------------------------------------------------------------+

Pipeline Execution
~~~~~~~~~~~~~~~~~~

The data pipeline can also be invoked from the command line.

.. code:: bash

   make load.database releases=<version> align=<boolean> kir=<boolean> limit=<int>

Retrieving logs, data and parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
See the reference section :ref:`makefilerefretrieve` for useful commands.

Developing Locally
------------------

.. note:: 
   This information is incomplete but will be updated soon.

Creating a Python Virtual Environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When developing locally, you will need to create an individual virtual
environment to run scripts in the ``jobs`` or ``functions`` directories,
since they require different dependencies.

.. code:: bash

   cd <specific job or function directory>
   python3 -m venv .venv
   source .venv/bin/activate
   pip install -U pip
   pip install -r requirements.txt

To use the virtual environment inside a Jupyter Notebook, first activate
the virtual environment, then create a kernel for it.

.. code:: bash

   # Install ipykernal
   pip install ipykernel python-dotenv

   # Add the kernel
   python3 -m ipykernel install --user --name=<environment name>

   # Remove the kernel
   jupyter kernelspec uninstall <environment name>
