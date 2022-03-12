[![Docker Image CI](https://github.com/databio/LOLAweb/actions/workflows/build.yml/badge.svg?branch=dev)](https://github.com/databio/LOLAweb/actions/workflows/build.yml)


# Shiny LOLAweb for Docker

The `ghcr.io/databio/lolaweb` container is based on the `ghcr.io/databio/shinybase` container, which you can find in its [github repository](https://github.com/databio/shinyBase) or in the [GitHub Container Registry](https://github.com/databio/shinyBase/pkgs/container/shinybase).

## `build` the container image yourself

1. Clone this repository
2. Build locally using Docker. Run this command from the same directory as the `Dockerfile`.

```docker build --no-cache -t lolaweb .```


## Or `pull` the container image:

```docker pull ghcr.io/databio/lolaweb```

The container image itself is hosted in the GitHub Container Registry: https://github.com/databio/LOLAweb/pkgs/container/lolaweb

## Container volumes and reference data

LOLAweb needs access to a few folders where it can store results or logs, or access necessary files like the database. To handle this, we've set up the app to look for two shell environment variables:

* **$LWREF**, for LOLAweb reference data, which may be read-only
* **$LWLOCAL**, where local results can be written.

To run the LOLAweb container (locally or on a server), you need to set these environment variables (for example, in your `.bashrc` or `.profile` file. These variables will be injected into the container when it is run.

For example, set these variables to two different paths if you like. Or if you keep all five subfolders together in the same path, set these variables to the same value.

```
# Example locations. Set to match your environment
export LREF='/home/subdir/lola'
export LWLOCAL='/var/things/loladata'
```


### LWREF

LOLAweb will look at the value in `$LWREF` for the reference data. This folder should have subfolders called `databases`, `universes`, and `examples`. In each of these subfolders are another layer of subfolders for genome assemblies. See the [LOLAweb documentation](https://github.com/databio/LOLAweb/tree/master/apps/LOLAweb) for downloading reference data.
 
### LWLOCAL

LOLAweb looks for `$LWLOCAL` to have two subfolders: `cache`, and `shinylogs`. This is where the app will write results and log files. If running LOLAweb on a server, be sure these directories are writeable by the Docker process.

## Run the LOLAweb container locally with reference data:

    docker run -d \
      -p 80:80 \
      -e LWREF=$LWREF \
      -e LWLOCAL=$LWLOCAL \
      --volume ${LWLOCAL}:${LWLOCAL} \
      --volume ${LWREF}:${LWREF} \
      --volume ${LWLOCAL}/shinylog:/var/log/shiny-server
      ghcr.io/databio/lolaweb

Open a browser to:
```
http://localhost/LOLAweb/apps/LOLAweb
```

## Running a `dev` container

You could also run the `dev` version of the container by pulling `ghcr.io/databio/lolaweb:dev`. This will retrieve the dev tagged image from the GitHub Container Registry. Just add `:dev` to the container name at the end of the `docker run` command above.


## Running multiple LOLAweb containers simultaneously with Docker Swarm:

For the typical use case of an individual user, a single running container will suffice. But if you need to set up an enterprise-level LOLAweb server that can handle concurrent users, we've also made that easy by using docker swarm. This is how we run the main LOLAweb servers, and you could do the same thing if you want your own local implementation. Docker Swarm is a technique for running multiple instances of the same container. [Read more](swarm/README.md) about how to set up your own swarm.
