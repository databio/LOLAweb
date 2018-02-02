FROM r-base:latest

MAINTAINER VP Nagraj "vpnagraj@virginia.edu"

# Install dependencies and Download and install shiny server
RUN apt-get update && apt-get install -y -t unstable \
    sudo \
    gdebi-core \
    pandoc \
    pandoc-citeproc \
    libsodium-dev \
    libssl-dev \
    libcurl4-gnutls-dev \
    libcairo2-dev/unstable \
    git \
    libxt-dev && \
    wget --no-verbose https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/VERSION -O "version.txt" && \
    VERSION=$(cat version.txt)  && \
    wget --no-verbose "https://s3.amazonaws.com/rstudio-shiny-server-os-build/ubuntu-12.04/x86_64/shiny-server-$VERSION-amd64.deb" -O ss-latest.deb && \
    gdebi -n ss-latest.deb && \
    rm -f version.txt ss-latest.deb && \
    R -e "install.packages(c('shiny', 'rmarkdown'), repos='https://cran.rstudio.com/')" && \
    cp -R /usr/local/lib/R/site-library/shiny/examples/* /srv/shiny-server/ && \
    rm -rf /var/lib/apt/lists/*

# install R packages
RUN R -e "install.packages(c('DT', 'ggplot2', 'shinyjs', 'devtools', 'shinyWidgets', 'sodium'), repos='http://cran.rstudio.com/')" -e "source('http://bioconductor.org/biocLite.R')" -e "biocLite(c('GenomicRanges', 'LOLA'))" -e "devtools::install_github('databio/simpleCache')"

# move conf files
COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY shiny-server.sh /usr/bin/shiny-server.sh

# get the code
WORKDIR /srv/shiny-server
RUN git clone https://github.com/databio/LOLAweb.git

# add dir for cache
RUN mkdir shinyLOLA/cache
RUN chown -R shiny:shiny shinyLOLA/cache

# run the server setup script
CMD ["/usr/bin/shiny-server.sh"]
