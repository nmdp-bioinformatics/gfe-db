Reference
=========

AWS Cloud Architecture
----------------------

.. image:: _static/img/arch-v220417.png
   :scale: 50%

``gfe-db`` architecture is organized into 3 layers.

#. `Base infrastructure <infrastructure_>`__ 
#. `Database <database_>`__ 
#. `Data Pipeline <datapipeline_>`__ 

This allows deployments to be decoupled using Makefiles. Common
configuration parameters are shared between resources using environment
variables, JSON files, AWS SSM Parameter Store and Secrets Manager.

.. _infrastructure:

Base infrastructure
~~~~~~~~~~~~~~~~~~~

The base infrastructure layer deploys a VPC, public subnet, S3 bucket,
Elastic IP and common SSM parameters and secrets for the other services
to use.

.. _database:

Database
~~~~~~~~

The database layer deploys an EC2 instance running the Neo4j Community
Edition into a public subnet. During initialization,
Cypher queries are run to create constraints and indexes, which help
speed up loading and ensure data integrity. Neo4j is ready to be
accessed through a browser once the instance has booted sucessfully.

.. _datapipeline:

Data Pipeline
~~~~~~~~~~~~~

The data pipeline layer automates integration of newly released IMGT/HLA
data into Neo4j using a scheduled Lambda which watches the source data
repository and invokes the build and load processes when it detects a
new IMGT/HLA version. The pipeline consists of a Step Functions state
machine which orchestrates two basic processes: build and load. The
build process employs a Batch job which produces an intermediate set of
CSV files. The load process leverages SSM Run Command to copy the CSV
files to the Neo4j server and execute Cypher statements directly on the
server (server-side loading). When loading the full dataset of 35,000+
alleles, the build step will generally take around 15 minutes, however
the load step can take an hour or more.
