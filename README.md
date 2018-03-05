# gfe-db
Graph database representing IPD-IMGT/HLA sequence data as GFE


## Building graph
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt

python bin/build_graph.py


## Docker
The easiest way to get the service running locally, is to pull an image containing the service from docker hub. Running the following command will pull the latest GFE service image from docker hub. The image on docker hub is built from the *Dockerfile* in the *docker* directory in the github repository. Every new commit to the *nmdp-bioinformatics/service-gfe-submission* repository triggers a new build of the docker image on docker hub.

```bash
docker pull nmdpbioinformatics/gfe-db:0.0.3260-hlai
docker run -d --name gfe-db -p 7474:7474 nmdpbioinformatics/gfe-db:0.0.3260-hlai
```
The *-d* flag runs the service in "detached-mode" in the background and *-p* specifies what ports to expose. Make sure the ports you expose are not already in use. If the docker container is successfuly executed then typing ``docker ps -a`` will show a new container labeled *service-gfe-submission* running. 

[Click here](https://hub.docker.com/r/nmdpbioinformatics/gfe-db/) for more information on the publically available docker image. 



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



