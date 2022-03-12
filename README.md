# LOLAweb

[![Docker Image CI](https://github.com/databio/LOLAweb/actions/workflows/build.yml/badge.svg?branch=dev)](https://github.com/databio/LOLAweb/actions/workflows/build.yml)

LOLAweb is a web server and interactive results viewer for enrichment of overlap between a query region set (a bed file) and a database of region sets. It provides an interactive result explorer to visualize the highest ranked enrichments from the database. You can access the web server at <http://lolaweb.databio.org>.

This repository contain the shiny [app source code](apps/LOLAweb/) and Docker implementation for LOLAweb. 

## Shiny app

LOLAweb is implemented as an interactive Shiny app. You can run this app locally by following the [instructions in the appfolder](apps/LOLAweb/).

## Docker

The `ghcr.io/databio/lolaweb` container is based on the `ghcr.io/databio/shinybase` container, which you can find in its [GitHub repository](https://github.com/databio/shinyBase) or in the [GitHub Container Registry](https://github.com/databio/shinyBase/pkgs/container/shinybase).

### `build` the container image yourself

1. Clone this repository
2. Build locally using Docker. Run this command from the same directory as the `Dockerfile`.

```docker build --no-cache -t lolaweb .```


### Or `pull` the container image:

```docker pull ghcr.io/databio/lolaweb```

The container image itself is hosted in the GitHub Container Registry: https://github.com/databio/LOLAweb/pkgs/container/lolaweb

### Container volumes and reference data

LOLAweb needs access to a few folders where it can store results or logs, or access necessary files like the database. To handle this, we've set up the app to look for two shell environment variables:

* **$LWREF**, for LOLAweb reference data, which may be read-only
* **$LWLOCAL**, where local results can be written.

To run the LOLAweb container (locally or on a server), you need to set these environment variables (for example, in your `.bashrc` or `.profile` file. These variables will be injected into the container when it is run.

For example, set these variables to two different paths if you like. Or if you keep all five subfolders together in the same path, set these variables to the same value.

```
# Example locations. Set to match your environment
export LWREF='/home/subdir/lola/'
export LWLOCAL='/var/things/loladata/'
```

LOLAweb will look at the value in `$LWREF` for the reference data. This folder should have subfolders called `databases`, `universes`, and `examples`. In each of these subfolders are another layer of subfolders for genome assemblies
 
LOLAweb looks for `$LWLOCAL` to have two subfolders: `cache`, and `shinylogs`. This is where the app will write results and log files. If running LOLAweb on a server, be sure these directories are writeable by the Docker process.

The following instructions demonstrat how to download and configure the LOLAweb data directories for a minimal example using `hg19` reference data:

```
## assign env vars for data path
## NOTE: must include trailing /
LWLOCAL="/path/to/local/data/"
LWREF="/path/to/reference/data/"

## change reference data dir
cd $LWREF

## create dir for databases
mkdir -p databases
## create examples and universe dir
## NOTE: these must include subdirs named corresponding to appropriate ref genome
mkdir -p examples/hg19
mkdir -p universes/hg19

## download example universe and user set
curl http://cloud.databio.org.s3.amazonaws.com/vignettes/lola_vignette_data_150505.tgz | tar xvz

## move example universe and user set files to hg19 dir
mv lola_vignette_data/activeDHS_universe.bed universes/hg19/.
mv lola_vignette_data/setB_100.bed examples/hg19/.

## clean up
rm -rf lola_vignette_data

## download databases
curl http://cloud.databio.org.s3.amazonaws.com/regiondb/LOLACoreCaches_170206.tgz | tar xvz
curl http://cloud.databio.org.s3.amazonaws.com/regiondb/LOLAExtCaches_170206.tgz | tar xvz

## move databases to appropriate spots
mv scratch/ns5bc/resources/regions/LOLACore databases/Core
mv scratch/ns5bc/resources/regions/LOLAExt databases/Extended

## clean up
rm -rf scratch

## change ot local data dir
cd $LWLOCAL

## create placeholder dirs for cache and shinylog
mkdir -p cache
mkdir -p shinylog
```

### Run the LOLAweb container locally with reference data:

```
## run the docker image
## NOTE: this run command uses image pulled from ghcr.io/databio/lolaweb
docker run -d \
  -p 80:80 \
  -e LWREF=$LWREF \
  -e LWLOCAL=$LWLOCAL \
  --volume ${LWLOCAL}:${LWLOCAL} \
  --volume ${LWREF}:${LWREF} \
  --volume ${LWLOCAL}/shinylog:/var/log/shiny-server \
  ghcr.io/databio/lolaweb
```

Open a browser to:
```
http://localhost/LOLAweb/apps/LOLAweb
```

### Running a `dev` container

You could also run the `dev` version of the container by pulling `ghcr.io/databio/lolaweb:dev`. This will retrieve the dev tagged image from the GitHub Container Registry. Just add `:dev` to the container name at the end of the `docker run` command above.

### Running multiple LOLAweb containers simultaneously with Docker Swarm:

For the typical use case of an individual user, a single running container will suffice. But if you need to set up an enterprise-level LOLAweb server that can handle concurrent users, we've also made that easy by using docker swarm. This is how we run the main LOLAweb servers, and you could do the same thing if you want your own local implementation. Docker Swarm is a technique for running multiple instances of the same container. [Read more](swarm/README.md) about how to set up your own swarm.

### Troubleshooting

The LOLAweb Docker implementation includes a mechanism to write Shiny Server logs to `$LWLOCAL/shinylog`. These log files may be useful when troubleshooting problems with running LOLAweb via Docker. They include errors with R processing as well as information as to whether the Shiny Server process was killed due to resource limitations (i.e., not enough RAM allocated to Docker daemon).

For additional support with the LOLAweb Docker implementation, please file a [GitHub issue](https://github.com/databio/LOLAweb/issues).

