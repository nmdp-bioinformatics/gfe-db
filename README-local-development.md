# Integrating repositories together
This file explains the integration of following repositories on local:

- [GFE-DB](https://github.com/nmdp-bioinformatics/gfe-db)
- [ACT-SERVICE](https://github.com/nmdp-bioinformatics/act-service)
- [ACT-CLIENT](https://github.com/nmdp-bioinformatics/act-client)
- [SEQ-ANN](https://github.com/nmdp-bioinformatics/seq-ann)
- [PY-ARD](https://github.com/nmdp-bioinformatics/py-ard)

### PY-ARD
The py-ard pip version is good enough to use, so this does not need update on local.

### GFE-DB
GFE-DB is the graph database representing IPD-IMGT/HLA sequence data as GFE. In order to run it locally:
- Pull the [source code](https://github.com/nmdp-bioinformatics/gfe-db).
- Go inside root folder and create a new empty directory named "output".
- Modify the line # 30 inside `bin/build.sh` and point it to your directories, in my case, I had to create a
an empty directory named "data" in project root and change this line to : `cp ./mod-imgt/C_gen.sth ./data/3360`
- You will be able to run following commands under "Build graph from source" successfully:
```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
cpanm install Bio::Perl
export IMGT="3360,3370"
export RELEASES=3360  # this value should be either 3360 or 3370 
export AN=True
bash bin/build.sh ./output
```
- Afterwards you will need to change the file `bin/load_graph.sh` and comment out line # 28 and uncomment
line # 27 to point it to your output directory.
- The file `bin/neo4j.conf.template` has a configuration template for gfe-db, this needs to be copied 
to `/opt/conf/neo4j.conf`. You can modify line # 50 in `bin/load_graph.sh` for this.
- Now you can run `bash bin/load_graph.sh  ./output`. You can see server running at `http://localhost:7474/`

### SEQ-ANN
The seqann allows users to annotate gene features in consensus sequences. In order to run it locally:
- Pull the [source code](https://github.com/nmdp-bioinformatics/seq-ann).
- Running the imgt_biosqldb from [DockerHub](https://hub.docker.com/r/nmdpbioinformatics/imgt_biosqldb/).
- Now try running the following code in python terminal:
```bash
from seqann import BioSeqAnn
from BioSQL import BioSeqDatabase
export BIOSQLHOST="localhost"
export BIOSQLPORT=3306
server = BioSeqDatabase.open_database(driver="pymysql", user="root",
                                      passwd="my-secret-pw", host="localhost",
                                      db="bioseqdb", port=3306)
seqann = BioSeqAnn(server=server)
ann = seqann.annotate(sequence, "HLA-A")
``` 
- seqann has some missing allele_lists at `seqann/data/allele_lists` . You need to download [Allelelist.3350.txt](https://raw.githubusercontent.com/ANHIG/IMGTHLA/3350/Allelelist.3350.txt), 
[Allelelist.3360.txt](https://raw.githubusercontent.com/ANHIG/IMGTHLA/3360/Allelelist.3360.txt) and 
[Allelelist.3370.txt](https://raw.githubusercontent.com/ANHIG/IMGTHLA/3370/Allelelist.3370.txt) to `seqann/data/allele_lists`.

This package is now ready to use with our integration process.

### ACT-SERVICE
The act-service is the API for annotation. In order to run it locally:
- Pull the [source code](https://github.com/nmdp-bioinformatics/act-service).
- Go inside repository root folder and inside `requirements.txt` remove the line # 38 "seqann".
- Now you should run following commands successfully:
```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
cd ../seq-ann/  # go to your local seq-ann package 
pip install .
cd ../act-service/  # go back to your act-service root folder.
```
- Go to file `swagger_server/controllers/type_align_controller.py` line # 24, update the `neo4j_url` to
`http://localhost:7474`.
- Now you need to update your port inside `main.py` (line # 19) to 80 and run the server again, you may 
need sudo permissions to do this

Your act-service API is good to go.

### ACT-CLIENT
The act-client is client for calling the ACT service. In order to run it locally:
- Pull the [source code](https://github.com/nmdp-bioinformatics/act-client).
- Now you should run following commands successfully:
```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```
- Go to `gfe_client/configuration.py` line # 50 and change it to `self.host = "http://localhost"`.
- Go to python terminal and run following: 
```bash
import gfe_client
api = gfe_client.TypeSeqApi()
seq = "AGAGACTCTCCCGAGGATTTCGTGTACCAGTTTAAGGCCATGTGCTACTTCACC"
response = api.typeseq_get(seq, imgthla_version="3.31.0", locus="HLA-A")
```
 - `response` will have the required result.
 
 
 This completes the integration process of above mentioned five repositories.