[![Docker pulls](https://img.shields.io/docker/pulls/somrc/shiny-lola.svg?style=plastic)](https://registry.hub.docker.com/v2/repositories/somrc/shiny-lola/) [![Docker Automated build](https://img.shields.io/docker/automated/somrc/shiny-lola.svg?style=plastic)](https://hub.docker.com/r/somrc/shiny-lola/)

# Shiny LOLAweb for Docker

## Build the container image yourself

1. Clone this repository
2. Build locally using Docker. Run this command from the same directory as the `Dockerfile`.

```docker build --no-cache -t shinylola .```


## Or Pull the container image:

```docker pull somrc/shiny-lola```


## Reference Data

See the [shinyLOLA documentation](https://github.com/databio/shinyLOLA) for retrieving reference data.


## Run the LOLA container locally with reference data:

    docker run -d \
      -p 80:80 \
      -v /home/lola/reference:/srv/shiny-server/shinyLOLA/reference \
      -v /home/lola/universes:/srv/shiny-server/shinyLOLA/universes \
      -v /home/lola/userSets:/srv/shiny-server/shinyLOLA/userSets \
      -v /home/lola/shinylog:/var/log/shiny-server \
      somrc/shiny-lola

Note that the fourth volume mounted is a log directory, allowing users to view logs generated from the Shiny application.

## Run LOLA on a server in a Docker Swarm:

Docker Swarm is a technique for running multiple instances of the same container.
[Read more](swarm/README.md) about how to set up your own swarm.
