Reference
=========

GFE Graph
------------

The GFE database (``gfe-db``) represents the relationships between GFEs,
features, sequences and other types of data. The new schema is centered
around the GFE node and makes the curation and database versioning of
WHO designations or WHO labels an optional annotation of GFEs.

.. image:: /_static/img/schema-alpha-v230305.png
   :scale: 50%
   :align: center

Breaking down a GFE
~~~~~~~~~~~~~~~~~~~

.. note::
   This section discusses the ``WHO`` and ``IMGT_HLA`` nodes which have 
   been deprecated.

The representation of a single GFE, for example corresponding to the
allele ``HLA-A*01:01:01:01`` can be understood from the graph.

The GFE node points to a WHO node as one of possibly many annotations.
With this schema it is possible to analyze GFEs that do not have an WHO
label associated with it.

To see the how a GFE expands to its constituent components, the
following query returns the corresponding features associated with the
GFE referred to by the WHO allele ``HLA-A*01:03:01:01``.

.. code::

   MATCH (:WHO {name:'HLA-A*01:03:01:01'})-[]-(:GFE)-[]-(f:Feature) 
   RETURN f.term, f.rank ORDER BY f.term, f.rank

Results:

+-----------------+--------+
| f.term          | f.rank |
+=================+========+
| EXON            | 1      |
+-----------------+--------+
| EXON            | 2      |
+-----------------+--------+
| EXON            | 3      |
+-----------------+--------+
| EXON            | 4      |
+-----------------+--------+
| EXON            | 5      |
+-----------------+--------+
| EXON            | 6      |
+-----------------+--------+
| EXON            | 7      |
+-----------------+--------+
| EXON            | 8      |
+-----------------+--------+
| FIVE_PRIME_UTR  | 1      |
+-----------------+--------+
| INTRON          | 1      |
+-----------------+--------+
| INTRON          | 2      |
+-----------------+--------+
| INTRON          | 3      |
+-----------------+--------+
| INTRON          | 4      |
+-----------------+--------+
| INTRON          | 5      |
+-----------------+--------+
| INTRON          | 6      |
+-----------------+--------+
| INTRON          | 7      |
+-----------------+--------+
| THREE_PRIME_UTR | 1      |
+-----------------+--------+

These features each have an accession number that is unique in the
context of the locus, term and rank combination and is a **permanent
reversible 1-to-1 mapping** between the sequence and the accession
number in that context.

Mapping is not permanent, reversible or 1-to-1 for these entity
properties: - WHO/IMGT_HLA names and sequences - WHO/IMGT_HLA names and
IMGT accession numbers - IMGT accession numbers and sequence

Here is a older example of a relationship between a WHO/IMGT_HLA allele
(``HLA-DRB1*11:17``) and the corresponding GFE.

.. image:: /_static/img/157B500B-9399-4B46-9E48-628F72C869C0.jpeg
   :scale: 50%
   :align: center

In this example, the GFE associated with this allele changed between
3.42.0 and 3.43.0

GFE nodes
~~~~~~~~~~~~~

Description
^^^^^^^^^^^

Each node represents a distinct GFE object. For example, a GFE with
``name="HLA-Aw2-1-1-1-1-4-1-1-1-2-1-1-1-1-1-1-4"`` corresponds to a full
sequence and also 17 features: - FIVE_PRIME_UTR - EXON (1-8) - INTRON
(1-7) - THREE_PRIME_UTR

Properties
^^^^^^^^^^

.. code:: json

   {
     "name": "HLA-Aw99-8-363-912-781-2901-128-581-151-324-198-9-316-80-508-43-30",
     "locus": "HLA-A"
   }


+--------------+--------------------------------------------------------------------+-----------+----------------------------------------+
| Property     | Example                                                            | Data Type | Description                            |
+==============+====================================================================+===========+========================================+
| ``name``     | HLA-Aw99-8-363-912-781-2901-128-581-151-324-198-9-316-80-508-43-30 | string    | GFE name                               |
+--------------+--------------------------------------------------------------------+-----------+----------------------------------------+
| ``locus``    | HLA-A                                                              | string    | Position of the gene on the chromosome |
+--------------+--------------------------------------------------------------------+-----------+----------------------------------------+

Feature nodes
~~~~~~~~~~~~~~~~~

.. _description-1:

Description
^^^^^^^^^^^

A feature is a tuple of: locus, term and rank. A locus is "anything in
`HUGO <https://www.genenames.org/>`__", and a term is "anything in
`sequence ontology <http://www.sequenceontology.org/>`__".

.. _properties-1:

Properties
^^^^^^^^^^

.. code:: json

   {
     "accession": 99,
     "locus": "HLA-A",
     "rank": 1,
     "sequence": "CAGGAGCAGAG...",
     "term": "FIVE_PRIME_UTR"
   }


+-------------+----------------+-----------+--------------------------------------------------------------------------------------------+
| Property    | Example        | Data Type | Description                                                                                |
+=============+================+===========+============================================================================================+
|``accession``| 2901           | string    | Relatively stable unique record identifier for a sequence                                  |
+-------------+----------------+-----------+--------------------------------------------------------------------------------------------+
|``locus``    | HLA-A          | string    | Position of the gene on the chromosome                                                     |
+-------------+----------------+-----------+--------------------------------------------------------------------------------------------+
|``rank``     | 7              | string    | Ordinal number describing the position of the Feature sequence on the allele               |
+-------------+----------------+-----------+--------------------------------------------------------------------------------------------+
|``sequence`` | CAGGAGCAGAG... | string    | Nucleotide sequence of the Feature                                                         |
+-------------+----------------+-----------+--------------------------------------------------------------------------------------------+
|``term``     | FIVE_PRIME_UTR | string    | Label describing the type of Feature; One of FIVE_PRIME_UTR, EXON, INTRON, THREE_PRIME_UTR |
+-------------+----------------+-----------+--------------------------------------------------------------------------------------------+

Sequence nodes
~~~~~~~~~~~~~~~~~~

.. _description-2:

Description
^^^^^^^^^^^

The nucleotide sequence corresponding to the GFE.

.. _properties-2:

Properties
^^^^^^^^^^

.. code:: json

   {
     "name": "HLA-Cw393-14-261-132-1610-454-45-532-107-272-205-3-264-71-398-4-621",
     "length": 3918,
     "locus": "HLA-C",
     "sequence": "TTATTTTGCTGGATGTAGTTTAATATTACCTGAGGTGAGGTAAGGTA..."
   }

+------------+---------------------------------------------------------------------+-----------+--------------------------------------------------------------------------+
| Property   | Example                                                             | Data Type | Description                                                              |
+============+=====================================================================+===========+==========================================================================+
|``name``    | HLA-Cw393-14-261-132-1610-454-45-532-107-272-205-3-264-71-398-4-621 | string    | Gene Feature Enumeration name                                            |
+------------+---------------------------------------------------------------------+-----------+--------------------------------------------------------------------------+
|``length``  | 3918                                                                | integer   | Length of nucleotide sequence                                            |
+------------+---------------------------------------------------------------------+-----------+--------------------------------------------------------------------------+
|``locus``   | HLA-C                                                               | string    | Position of the gene on the chromosome                                   |
+------------+---------------------------------------------------------------------+-----------+--------------------------------------------------------------------------+
|``sequence``| TTATTTTGCTGGATGTAGTTTAATATTACCTGAGGTGAGGTAAGGTA...                  | string    | Full nucleotide sequence                                                 |
+------------+---------------------------------------------------------------------+-----------+--------------------------------------------------------------------------+

IPD_Allele nodes
~~~~~~~~~~~~~~~~

.. _description-3a:

Description
^^^^^^^^^^^

Represents the IPD-IMGT/HLA designation for the allele.

.. _properties-3a:

Properties
^^^^^^^^^^

.. code:: json

   {
     "gene": "HLA-A",
     "lg": "HLA-A*01:242g",
     "name": "HLA-A*01:242" 
   }


+----------+-------------------+-----------+------------------------+
| Property | Example           | Data Type | Description            |
+==========+===================+===========+========================+
|``gene``  | HLA-A             |string     | Name of gene           |
+----------+-------------------+-----------+------------------------+
|``lg``    | HLA-A*01:242g     |string     | Locus group            |
+----------+-------------------+-----------+------------------------+
|``name``  | HLA-A*01:242      |string     | Name of allele         |
+----------+-------------------+-----------+------------------------+

IPD_Accession nodes
~~~~~~~~~~~~~~~~~~~

.. _description-3b:

Description
^^^^^^^^^^^

Represents the IPD-IMGT/HLA accession number for the allele.

.. _properties-3b:

Properties
^^^^^^^^^^

.. code:: json

   {
     "name": "HLA17766"
   }


+----------+-------------------+-----------+------------------------+
| Property | Example           | Data Type | Description            |
+==========+===================+===========+========================+
|``name``  | HLA17766          |string     | Name of allele         |
+----------+-------------------+-----------+------------------------+


Submitter nodes
~~~~~~~~~~~~~~~~~~~

.. _description-4:

Description
^^^^^^^^^^^

Describes the submitter of a GFE node and their contact and organization information.

.. _properties-4:

Properties
^^^^^^^^^^

.. code:: json

   {
     "email": "<email>",
     "institution": "<institution name>",
     "name": "<name>",
     "url": "<url>"
   }


+---------------+------------------------------------+-----------+---------------------------+
| Property      | Example                            | Data Type | Description               |
+===============+====================================+===========+===========================+
|``email``      | <email>                            | string    | Submitter's email         |
+---------------+------------------------------------+-----------+---------------------------+
|``institution``| IPD                                | integer   | Submitter's institution   |
+---------------+------------------------------------+-----------+---------------------------+
|``name``       | IPD-IMGT                           | string    | Submitter's full name     |
+---------------+------------------------------------+-----------+---------------------------+
|``url``        | https://www.ebi.ac.uk/ipd/imgt/hla/| string    | Submitter's website       |
+---------------+------------------------------------+-----------+---------------------------+


HAS_FEATURE edges
~~~~~~~~~~~~~~~~~~~~

.. _description-5:

Description
^^^^^^^^^^^

Links a GFE node to a Feature node.

.. _properties-5:

Properties
^^^^^^^^^^

No properties.

HAS_SEQUENCE edges
~~~~~~~~~~~~~~~~~~~~~

.. _description-6:

Description
^^^^^^^^^^^

Links a GFE node to the full Sequence node.

.. _properties-6:

Properties
^^^^^^^^^^

No properties.

HAS_IPD_Allele edges
~~~~~~~~~~~~~~~~~~~~

.. _description-7a:

Description
^^^^^^^^^^^

Links a GFE node to the IPD_Allele node.

.. _properties-7a:

Properties
^^^^^^^^^^

.. code:: json

   {
     "releases": [3500, 3510]
   }


+------------+-------------------+----------------+----------------------------------------------+
| Property   | Example           | Data Type      | Description                                  |
+============+===================+================+==============================================+
|``releases``| [..., 3500, 3510] | array[integer] | Release versions containing the relationship |
+------------+-------------------+----------------+----------------------------------------------+

HAS_IPD_ALLELE edges
~~~~~~~~~~~~~~~~~~~~

.. _description-7b:

Description
^^^^^^^^^^^

Links an IPD_Allele node to the IPD_Accession node.

.. _properties-7b:

Properties
^^^^^^^^^^

.. code:: json

   {
     "releases": [3500, 3510]
   }


+------------+-------------------+----------------+----------------------------------------------+
| Property   | Example           | Data Type      | Description                                  |
+============+===================+================+==============================================+
|``releases``| [..., 3500, 3510] | array[integer] | Release versions containing the relationship |
+------------+-------------------+----------------+----------------------------------------------+


SUBMITTED edges
~~~~~~~~~~~~~~~~~~

.. _description-8:

Description
^^^^^^^^^^^

Links the Submitter node to the GFE node.

.. _properties-8:

Properties
^^^^^^^^^^

.. code:: json

   {
     "submit_date": "2022-02-17"
   }


+---------------+------------+-----------------+--------------------+
| Property      | Example    | Data Type       | Description        |
+===============+============+=================+====================+
|``submit_date``| 2022-02-17 | datetime string | Date of submission |
+---------------+------------+-----------------+--------------------+

Service Configurations
----------------------

Configuring is managed using JSON files, SSM Parameter Store, Secrets
Manager, and shell variables. To deploy changes in these files, run the
command.

.. code:: bash

   make config.deploy

Graph Database
~~~~~~~~~~~~~~

Neo4j
^^^^^

Custom configuration settings for Neo4j are contained in
``neo4j.template``. This file is copied into ``/etc/neo4j`` during boot
or manually. When Neo4j is restarted it will use the settings in
``neo4j.template`` to overwrite ``neo4j.conf``. More information can be
found in the documentation here at `Neo4j Cloud Virtual Machines
<https://neo4j.com/developer/neo4j-cloud-vms/>`_.

.. important::
   Neo4j no longer supports the Community Edition of their AMI for EC2.
   The next release of ``gfe-db`` will use the Bitnami Neo4j AMI which 
   will change this information.

Shell Scripts
^^^^^^^^^^^^^

Bash scripts are used for automating Neo4j configuration, loading and
backup. These are stored in S3 and executed on the database instance using 
SSM Run Command. These are found in ``gfe-db/gfe-db/database/scripts/``.

To update shell scripts on the Neo4j instance, run the following commands in sequence.

.. code:: bash

    # sync the scripts to S3
    make config.deploy

    # sync the scripts from S3 to the instance
    make database.sync-scripts

Cypher Scripts
^^^^^^^^^^^^^^

Cypher scripts manage node constraints & indexes and load the data.
These are found in ``gfe-db/gfe-db/database/neo4j/cypher/``.

.. _datapipelineconfig:

Data Pipeline
~~~~~~~~~~~~~

Input Parameters
^^^^^^^^^^^^^^^^

The ``pipeline-input.json`` is stored in S3 and contains the default
configuration used for automated updates.

.. code:: json

   // pipeline-input.json
   {
     "align": "False",
     "kir": "False",
     "mem_profile": "False",
     "limit": ""
   }

IMGT/HLA Release Versions State
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The application’s state tracks which releases have been processed and
added to the database. This file tracks the releases which have already
been processed. If the ``gfe-db-invoke-pipeline`` function detects a
valid release branch in the source data repository that is not in the
``releases`` array, it will start the pipeline for this release. Once
the update is finished, the processed release is appended to the array.

.. code:: json

   // IMGTHLA-repository-state.json
   {
     "timestamp": "2021-12-09 02:36:59",
     "repository_url": "https://github.com/ANHIG/IMGTHLA",
     "releases": [
       "3100",
       // ...,
       "3470"
     ]
   }

+----------------+----------------------------------+------------------+------------------------------------------------------------------+
| Variable       | Example Value                    | Type             | Description                                                      |
+================+==================================+==================+==================================================================+
| repository_url | https://github.com/ANHIG/IMGTHLA | string           | The repository the trigger is watching                           |
+----------------+----------------------------------+------------------+------------------------------------------------------------------+
| releases       | ["3100", ..., "3470"]            | array of strings | List of available releases. Any release added to the repository  |
|                |                                  |                  | that is not in this list will trigger the pipeline build.        |
+----------------+----------------------------------+------------------+------------------------------------------------------------------+

Logging
~~~~~~~
Logs for EC2, Lambda and Batch are collected by CloudWatch Logs. 

.. _makefileref:

Makefile Command Reference
--------------------------

To see a list of possible commands using Make, run ``make`` on the
command line.

Deploy to AWS
~~~~~~~~~~~~~

Deploy all CloudFormation based services:

.. code:: bash

   make deploy

Deploy specific stacks.

.. code:: bash

   make infrastructure.deploy

.. code:: bash

   make database.deploy

.. code:: bash

   make pipeline.deploy

Load releases
~~~~~~~~~~~~~

Run the StepFunctions State Machine to load Neo4j:

.. code:: bash

   make database.load.run releases=<version> align=<boolean> kir=<boolean> limit=<int>

.. _makefilerefretrieve:

Manage Services
~~~~~~~~~~~~~~~

.. code:: bash

   make database.start

.. code:: bash

   make database.stop

.. code:: bash

   make database.status

.. code:: bash

   make database.reboot

Backup & restore
~~~~~~~~~~~~~~~~

Backup the database to S3:

.. code:: bash

   make database.backup

List the available backups:

.. code:: bash

   make database.backup.list

Restore from a backup:

.. code:: bash

   make database.restore from_date=<YYYY/MM/DD/HH>


Sync and retrieve logs, data and configuration values
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Deploy config files and scripts to S3:

.. code:: bash

   make config.deploy

Sync the Neo4j database scripts to S3:

.. code:: bash

   make database.sync-scripts

Download CSV data from builds to S3 to ``./data``:

.. code:: bash

   make get.data

Download logs from EC2 to ``./logs``:

.. code:: bash

   make get.logs

Display the Neo4j Browser endpoint URL:

.. code:: bash

   make database.get.endpoint

Fetch the database credentials:

.. code:: bash

   make database.get.credentials

Fetch the instance ID:

.. code:: bash

   make database.get.instance-id

Tear down infrastructure
~~~~~~~~~~~~~~~~~~~~~~~~

Delete all CloudFormation based services and data:

.. code:: bash

   make delete

Delete specific stacks (may cause issues due to missing or stale dependencies):

.. code:: bash

   make infrastructure.delete

.. code:: bash

   make database.delete

.. code:: bash

   make pipeline.delete

Build documentation
~~~~~~~~~~~~~~~~~~~

.. code:: bash

   make docs.build

View the docs in a web browser:

.. code:: bash

   make docs.url

Backup & Restore
----------------

Backups
~~~~~~~

Backups are orchestrated by Systems Manager. To create a backup, run the command.

.. code:: bash

   make database.backup

This will create a backup of the Neo4j database and store it in S3 under the path
``s3://<data bucket name>/backups/neo4j/YYYY/MM/DD/HH/gfedb.zip``.

To see a list of available backup dates that can be restored, run the command.

.. code:: bash

   make database.backup.list

Restore
~~~~~~~

To restore from a backup, run the command and pass the date of the backup you wish to restore
using the format YYYY/MM/DD/HH.

.. code:: bash

   make database.restore from_date=<YYYY/MM/DD/HH>