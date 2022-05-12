Reference
=========

GFE Graph
------------

The GFE database (``gfe-db``) represents the relationships between GFEs,
features, sequences and other types of data. The new schema is centered
around the GFE node and makes the curation and database versioning of
WHO designations or WHO labels an optional annotation of GFEs.

.. image:: /_static/img/schema-alpha-v220511.png
   :scale: 50%

GFE nodes
~~~~~~~~~~~~~

Description
^^^^^^^^^^^

Each node represents a distinct GFE object. For example, a GFE with
``gfe_name="HLA-Aw2-1-1-1-1-4-1-1-1-2-1-1-1-1-1-1-4"`` corresponds to a full
sequence and also 17 features: - FIVE_PRIME_UTR - EXON (1-8) - INTRON
(1-7) - THREE_PRIME_UTR

Properties
^^^^^^^^^^

.. code:: json

   {
     "gfe_name": "HLA-Aw99-8-363-912-781-2901-128-581-151-324-198-9-316-80-508-43-30",
     "locus": "HLA-A"
   }


+--------------+--------------------------------------------------------------------+-----------+----------------------------------------+
| Property     | Example                                                            | Data Type | Description                            |
+==============+====================================================================+===========+========================================+
| ``gfe_name`` | HLA-Aw99-8-363-912-781-2901-128-581-151-324-198-9-316-80-508-43-30 | string    | GFE name                               |
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
     "gfe_name": "HLA-Cw393-14-261-132-1610-454-45-532-107-272-205-3-264-71-398-4-621",
     "length": 3918,
     "locus": "HLA-C",
     "seq_id": 27670532806245477286153332635897,
     "sequence": "TTATTTTGCTGGATGTAGTTTAATATTACCTGAGGTGAGGTAAGGTA..."
   }

+------------+---------------------------------------------------------------------+-----------+--------------------------------------------------------------------------+
| Property   | Example                                                             | Data Type | Description                                                              |
+============+=====================================================================+===========+==========================================================================+
|``gfe_name``| HLA-Cw393-14-261-132-1610-454-45-532-107-272-205-3-264-71-398-4-621 | string    | Gene Feature Enumeration name                                            |
+------------+---------------------------------------------------------------------+-----------+--------------------------------------------------------------------------+
|``length``  | 3918                                                                | integer   | Length of nucleotide sequence                                            |
+------------+---------------------------------------------------------------------+-----------+--------------------------------------------------------------------------+
|``locus``   | HLA-C                                                               | string    | Position of the gene on the chromosome                                   |
+------------+---------------------------------------------------------------------+-----------+--------------------------------------------------------------------------+
|``seq_id``  | 27670532806245477286153332635897                                    | integer   | Compressed UUID based on MD5 hash of sequence (used for faster indexing) |
+------------+---------------------------------------------------------------------+-----------+--------------------------------------------------------------------------+
|``sequence``| TTATTTTGCTGGATGTAGTTTAATATTACCTGAGGTGAGGTAAGGTA...                  | string    | Full nucleotide sequence                                                 |
+------------+---------------------------------------------------------------------+-----------+--------------------------------------------------------------------------+

IPD_Allele nodes
~~~~~~~~~~~~~~~~

.. note::
   ``IPD_Allele`` and ``IPD_ACC`` nodes replace the previous ``WHO`` nodes. Documentation
   is in progress.

.. _description-3a:

Description
^^^^^^^^^^^

*Documentation in progress*

.. _properties-3a:

Properties
^^^^^^^^^^

.. code:: json

   {
     // Documentation in progress
   }

+----------+-------------------+-----------+------------------------+
| Property | Example           | Data Type | Description            |
+==========+===================+===========+========================+
|          |                   |           |                        |
+----------+-------------------+-----------+------------------------+

IPD_ACC nodes
~~~~~~~~~~~~~

.. note::
   ``IPD_Allele`` and ``IPD_ACC`` nodes replace the previous ``WHO`` nodes. Documentation
   is in progress.

.. _description-3b:

Description
^^^^^^^^^^^

*Documentation in progress*

.. _properties-3b:

Properties
^^^^^^^^^^

.. code:: json

   {
     // Documentation in progress
   }

+----------+-------------------+-----------+------------------------+
| Property | Example           | Data Type | Description            |
+==========+===================+===========+========================+
|          |                   |           |                        |
+----------+-------------------+-----------+------------------------+

Submitter nodes
~~~~~~~~~~~~~~~~~~~

.. _description-4:

Description
^^^^^^^^^^^

Describes the submitter of a GFE node.

.. _properties-4:

Properties
^^^^^^^^^^

.. code:: json

   {
     "email": "<email>",
     "institution": "<institution name>",
     "name": "<name>"
   }

+---------------+----------------------+-----------+-------------------------+
| Property      | Example              | Data Type | Description             |
+===============+======================+===========+=========================+
|``email``      | @cibmtr.org          | string    | Submitter's email       |
+---------------+----------------------+-----------+-------------------------+
|``institution``| CIBMTR               | integer   | Submitter's institution |
+---------------+----------------------+-----------+-------------------------+
|``name``       | first name last name | string    | Submitter's full name   |
+---------------+----------------------+-----------+-------------------------+

HAS_FEATURE edges
~~~~~~~~~~~~~~~~~~~~

.. _description-5:

Description
^^^^^^^^^^^

Links a GFE node to a Feature node.

.. _properties-5:

Properties
^^^^^^^^^^

.. code:: json

   {
     // No properties
   }

HAS_SEQUENCE edges
~~~~~~~~~~~~~~~~~~~~~

.. _description-6:

Description
^^^^^^^^^^^

Links a GFE node to the full Sequence node.

.. _properties-6:

Properties
^^^^^^^^^^

.. code:: json

   {
     // No properties
   }

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
     "releases": [3470, 3460]
   }

+------------+--------------+----------------+----------------------------------------------+
| Property   | Example      | Data Type      | Description                                  |
+============+==============+================+==============================================+
|``releases``| [3470, 3460] | array[integer] | Release versions containing the relationship |
+------------+--------------+----------------+----------------------------------------------+

HAS_IPD_Allele edges
~~~~~~~~~~~~~~~~~~~~

.. _description-7b:

Description
^^^^^^^^^^^

Links an IPD_Allele node to the IPD_ACC node.

.. _properties-7b:

Properties
^^^^^^^^^^

.. code:: json

   {
     "releases": [3470, 3460]
   }

+------------+--------------+----------------+----------------------------------------------+
| Property   | Example      | Data Type      | Description                                  |
+============+==============+================+==============================================+
|``releases``|  3470        |        integer | Release versions containing the relationship |
+------------+--------------+----------------+----------------------------------------------+

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

Breaking down a GFE
~~~~~~~~~~~~~~~~~~~

.. deprecated:: 0.1
   This section discusses the ``WHO`` node which is has been deprecated. 
   We hope to update this as quickly as possible. 

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

In this example, the GFE associated with this allele changed between
3.42.0 and 3.43.0

.. _architecture:

AWS Cloud Architecture
----------------------

.. image:: /_static/img/arch-v220417.png
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