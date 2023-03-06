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

#. Purchase or designate a domain in Route53 and create a hosted zone
#. Acquire a subscription for the Bitnami Neo4j AMI through AWS Marketplace
#. Install `prerequisites <prerequisites_>`__ 
#. Set `environment variables <environment_>`__ 
#. Check the config JSONs (parameters and state) and edit the values as desired 
#. Run ``make deploy`` to deploy the stacks to AWS 
#. Run ``database.load.run releases=<version>`` to load Neo4j, or ``make database.load.run releases=<version> limit=<limit>``
#. Run ``make database.get-credentials`` to get the username and password for Neo4j
#. Run ``make database.get-url`` to get the URL for Neo4j and navigate to the Neo4j browser at the subdomain and host domain, for example ``https://gfe-db.cloudftl.com:7473/browser/``

.. _prerequisites:

Prerequisites
~~~~~~~~~~~~~

Please refer to the respective documentation for specific installation
instructions.

* Route53 domain and hosted zone
* Bitnami Neo4j AMI subscription and AMI ID
* GNU Make 3.81
* coreutils (optional but recommended)
* AWS CLI
* SAM CLI
* Docker
* jq
* Python 3.9 (if developing locally)

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
    GITHUB_PERSONAL_ACCESS_TOKEN=<secret>
    HOST_DOMAIN=<Fully qualified domain name>
    SUBDOMAIN=<subdomain>
    ADMIN_EMAIL=<email>
    APOC_VERSION=4.4.0.3
    GDS_VERSION=2.0.1
    NEO4J_AMI_ID=ami-04aa5da301f99bf58 # Bitnami Neo4j image available through AWS Marketplace

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

Run the command to deploy the application to AWS.

.. code:: bash

   make deploy

It is also possible to deploy or update individual services on at a time.

.. code:: bash

   # Deploy/update only the database service
   make database.deploy

   # Deploy/update only the pipeline service
   make pipeline.deploy

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

   make config.deploy

Clean Up
~~~~~~~~

To tear down resources run the command. You will need to manually delete
the data in the S3 bucket first.

.. code:: bash

   make delete

Use the following commands to tear down individual services.

.. code:: bash

   # Delete only the database service
   make database.delete

   # Delete only the pipeline service
   make pipeline.delete

.. warning::
   Deleting and re-deploying a layer may cause the parameters shared by
   other layers to go out of date. To avoid this the recommendation is to deploy
   sequentially and teardown in reverse sequence. For example, tearing down and
   re-deploying the database stack may affect parameters shared with the pipeline stack, 
   causing the pipeline stack to fail. The solution would be to also be 
   tear down and re-deploy the pipeline stack.

Loading Releases
----------------
Once deployed, the pipeline will check the source repository for new releases. If a new release is found the pipeline 
will trigger an execution to build and load the data into Neo4j.

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
| LIMIT       | null          | string | Number of alleles to build. Leave blank ("") to build all alleles. |
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

   database.load.run releases=<version> align=<boolean> kir=<boolean> limit=<int>

Retrieving logs, data and parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
See the reference section :ref:`makefilerefretrieve` for useful commands.

Developing Locally
------------------

When developing locally, you will need to create an individual virtual
environment to run scripts in the ``jobs`` or ``functions`` directories,
since they require different dependencies.

.. code:: bash

   cd <specific job or function directory>
   python3 -m venv .venv-dev
   source .venv-dev/bin/activate
   pip install -U pip
   pip install -r requirements-dev.txt

To use the virtual environment inside a Jupyter Notebook, first activate
the virtual environment, then create a kernel for it.

.. code:: bash

   # Install ipykernel
   pip install ipykernel python-dotenv

   # Add the kernel
   python3 -m ipykernel install --user --name=<environment name>

   # Remove the kernel
   jupyter kernelspec uninstall <environment name>

The build job requires a specific environment which can be created wth the following steps.
1. Add the following variables to `.env` inside the `build/` directory in the pipeline layer.

.. code:: bash

    REGION=<AWS region>
    GFE_BUCKET=<name of the deployed S3 bucket>
    FAILED_ALLELES_QUEUE=<name of the deployed Failed Alleles queue>
    ALIGN=False
    KIR=False
    MEM_PROFILE=True
    RELEASES=<release>
    LIMIT=<int>

2. Source the variables to the environment.

.. code:: bash

   # Export .env file variables
   set -a && source .env.dev && set +a

   # Check that the variables were set
   env