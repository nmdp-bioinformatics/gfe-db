# gfe-db
Graph database representing IPD-IMGT/HLA sequence data as GFE


## Docker
The easiest way to get the service running locally, is to pull an image containing the service from Docker hub. You can build the image locally and increase the number of IMGT/HLA releases loaded and include KIR or pull the pre-built image from Docker hub. The environment variable **RELEASES** specifies how many IMGT/HLA release you want to be loaded into the graph. The deafult number of releases is one. The environment variable **KIR** specifies if you want KIR to be loaded or not. By default KIR data is not loaded into the graph.

#### Building and running locally
```bash
git clone https://github.com/nmdp-bioinformatics/gfe-db
cd gfe-db
# Builds the graph with 3 IMGT/HLA DB versions and also adds KIR data
docker build -t gfe-db --build-arg IMGT="3310,3300" --build-arg AN=True .
docker run -d --name gfe-db -p 7474:7474 gfe-db
```

#### Pulling from Docker hub
```bash
# ** Image on Docker hub has only 1 IMGT release and KIR data is not loaded **
docker pull nmdpbioinformatics/gfe-db
docker run -d --name gfe-db -p 7474:7474 nmdpbioinformatics/gfe-db
```
The *-d* flag runs the service in "detached-mode" in the background and *-p* specifies what ports to expose. Make sure the ports you expose are not already in use. If the docker container is successfuly executed then typing ``docker ps -a`` will show a new container labeled **gfe-db** running. 

[Click here](https://hub.docker.com/r/nmdpbioinformatics/gfe-db/) for more information on the publically available docker image. 

## Building graph from source
The environment variable RELEASES specifies how many IMGT/HLA release you want to be loaded into the graph.

```
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
cpanm install Bio::Perl
export IMGT="3360,3370"
export AN=True
sh bin/build.sh /output/directory
sh bin/load_graph.sh  /output/directory
```


### Related Links

 * [hub.docker.com/r/nmdpbioinformatics/service-gfe-submission](https://hub.docker.com/r/nmdpbioinformatics/service-gfe-submission)
 * [service-gfe-submission.readthedocs.io](https://service-gfe-submission.readthedocs.io/en/latest/index.html)
 * [github.com/nmdp-bioinformatics/service-feature](https://github.com/nmdp-bioinformatics/service-feature)
 * [github.com/nmdp-bioinformatics/HSA](https://github.com/nmdp-bioinformatics/HSA)
 * [bioinformatics.bethematchclinical.org](https://bioinformatics.bethematchclinical.org)
 * [feature.nmdp-bioinformatics.org](https://feature.nmdp-bioinformatics.org)
 * [gfe.b12x.org](http://gfe.b12x.org)

<p align="center">
  <img src="https://bethematch.org/content/site/images/btm_logo.png">
</p>



