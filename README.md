# gfe-db
Graph database representing IPD-IMGT/HLA sequence data as GFE.

## Running the GFE database in Neo4j 4.2 using Docker
This README outlines the steps for building and running a development version of `gfe-db` in a local Docker container. Docker will deploy an instance of Neo4j 4.2 including the [APOC](https://neo4j.com/labs/apoc/4.1/) and [Graph Data Science](https://neo4j.com/docs/graph-data-science/current/) plugins. GFE data is stored in the `data/csv/` directory which is mounted as an external volume within the container when run. This keeps the data outside the container so that it can be updated easily.

## Development Constraints
This pipeline is under development.
* The graph schema is being redesigned.
* KIR data is not yet included.
* Releases other than IPD-IMGT/HLA Release 3.36.0 are not yet included (coming next).
* Loading the full GFE dataset of 20,000+ alleles might fail, or take an excessive amount of time. (Load scripts are in the processed of being optimized for this.)

Please feel free to open issues regarding specific bugs and feature requests.

## Project Files
```bash
.
├── bin                             # Executable scripts
│   ├── __init__.py
│   ├── build.sh                    # Entrypoint for build step
│   ├── build_gfedb.py              # Generates CSVs for Neo4j graph
│   └── get_alignments.sh           # Alignments are included by default            
├── (data)                          # Created during build step
│   ├── 3360                        # IMGT/HLA release 3.36.0
│   ├── csv                         # CSVs loaded into Neo4j
│   │   ├── all_alignments.3360.csv
│   │   ├── all_cds.3360.csv
│   │   ├── all_features.3360.csv
│   │   ├── all_groups.3360.csv
│   │   └── gfe_sequences.3360.csv
│   └── hla.3360.dat                # Allele data
├── neo4j                           # Neo4j load scripts
│   ├── bulk_load.cyp               # Under development
│   └── load.cyp                    # Merges new nodes with existing
├── notebooks                       # Development jupyter notebooks
├── .dockerignore                   # Files for Docker to ignore
├── .gitignore                      # Files for git to ignore
├── Dockerfile                      # Docker image for Neo4j 4.2
├── LICENSE
├── README-neo4j.md                 # Instructions for this workflow
├── README.md                       # gfe-db README
└── requirements.txt                # Python dependencies
```
## 1. Getting Started
Clone the repo.
```bash
git clone https://github.com/abk7777/gfe-db.git
```
Create a virtual environment and activate.
```bash
# Create .venv
python3 -m venv .venv

# Activate
source .venv/bin/activate
```
Install the requirements.
```bash
pip install -r requirements.txt
```

## 2. Build the GFE dataset
Run this script to generate a set CSV files of GFE data in the `data/csv/` directory. It is recommended to limit the number of alleles to avoid excessive build and load times.
```bash
# Limit the build to 100 alleles (recommended)
bash bin/build.sh 100

# Build complete database (not recommended)
bash bin/build.sh
```

## 3. Build Neo4j Docker image
Build the Docker image as defined in the Dockerfile. See [Configuring Neo4j in Dockerfile](#Configuring-Neo4j-in-Dockerfile) for important configuration settings.
```
docker build --tag gfe-db .
```

## 4. Start Neo4j graph database
Run the container to start Neo4j in Docker.
```
# Run container to start Neo4j
docker run -d --name gfe -v "$(pwd)"/data/csv/:/var/lib/neo4j/import \
    -p 7474:7474 -p 7473:7473 -p 7687:7687 gfe-db
```
If desired, access the container logs during startup. This will indicate when Neo4j is ready.
```bash
docker logs -f gfe
```
Stop and restart when needed.
```bash
# Stop container
docker stop gfe

# Start container
docker start gfe
```
## 5. Load the GFE data
Once the container is running and the Neo4j server is up, the data can be loaded using the Cypher script.
```
cat neo4j/load.cyp | docker exec --interactive gfe cypher-shell -u neo4j -p gfedb
```
*Note: This step is not yet optimized for the full dataset, so proceed with caution. For local development on the GFE graph, it is recommended to specify a limited number of alleles during the build step.*
## 6. Access Neo4j
Neo4j can be accessed through web browser at [http://localhost:7474/browser/](http://localhost:7474/browser/) and the data can be queried using Cypher.

To view the schema, run this command.
```cypher
CALL db.schema.visualization;
```

## 7. Clean up Docker
Delete the Docker container.
```bash
docker stop gfe
docker rm gfe
```
Delete the Docker image.
```bash
# List the images and get the IMAGE IDs for gfe-db and neo4j
docker image ls

# Remove the images for gfe-db:latest and neo4j:4.2 by the IMAGE ID
docker image rm <IMAGE ID> <IMAGE ID>
```

The fastest way to remove all Docker images, containers and volumes is the `prune` method. *Use with caution because this will delete all Docker images, containers and their data on your machine.*
```bash
# Use with caution
docker system prune --volumes -a
```
# Configuring Neo4j in Dockerfile
Configuration settings for Neo4j are passed through environment variables in the Dockerfile.
## Username & Password
The username and password is set as follows:
```Dockerfile
# Dockerfile
ENV NEO4J_AUTH=neo4j/gfedb
```
## Memory Management
Optimal memory for Neo4j depends on available RAM. Loading and querying a larger dataset will require more memory allocated. For more information on memory management in Neo4j, see the [Neo4j Operations Manual](https://neo4j.com/docs/operations-manual/current/performance/memory-configuration/).
```Dockerfile
# Dockerfile
ENV NEO4J_dbms_memory_heap_initial__size=2G \
    NEO4J_dbms_memory_heap_max__size=2G

# Rebuild the image after updating these
```

# Related Links

 * [hub.docker.com/r/nmdpbioinformatics/service-gfe-submission](https://hub.docker.com/r/nmdpbioinformatics/service-gfe-submission)
 * [service-gfe-submission.readthedocs.io](https://service-gfe-submission.readthedocs.io/en/latest/index.html)
 * [github.com/nmdp-bioinformatics/service-feature](https://github.com/nmdp-bioinformatics/service-feature)
 * [github.com/nmdp-bioinformatics/HSA](https://github.com/nmdp-bioinformatics/HSA)
 * [bioinformatics.bethematchclinical.org](https://bioinformatics.bethematchclinical.org)
 * [feature.nmdp-bioinformatics.org](https://feature.nmdp-bioinformatics.org)
 * [gfe.b12x.org](http://gfe.b12x.org)

<br>

<p align="center">
  <img src="https://bethematch.org/content/site/images/btm_logo.png">
</p>