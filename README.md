# LOLAweb
LOLAweb is a web server and interactive results viewer for enrichment of overlap between a query region set (a bed file) and a database of region sets. It provides an interactive result explorer to visualize the highest ranked enrichments from the database. You can access the web server at <http://lolaweb.databio.org>.

This repository contains two components: 1) the shiny [app source code](apps/LOLAweb/) and 2) [Docker implementation](docker/) for LOLAweb. 

## Shiny app

LOLAweb is implemented as an interactive shiny app. You can run this app locally by following the [instructions in the appfolder](apps/LOLAweb/).

## Dockerfile

We have also produced a Dockerfile and host an image on dockerhub (`databio/lolaweb`). You can read how to run LOLAweb locally in a docker container by reading the [instructions in the docker folder](docker/). 
