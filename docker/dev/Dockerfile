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
    libxml2-dev \
    libmariadb-dev \
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
RUN R -e "install.packages(c('DT', 'ggplot2', 'shinyjs', 'devtools', 'shinyWidgets', 'sodium', 'shinyBS'), repos='http://cran.rstudio.com/')" -e "source('http://bioconductor.org/biocLite.R')" -e "biocLite(c('GenomicRanges', 'LOLA', 'IRanges', 'ensembldb', 'BSgenome', 'BSgenome.Hsapiens.UCSC.hg38.masked', 'BSgenome.Hsapiens.UCSC.hg19.masked', 'BSgenome.Hsapiens.UCSC.mm10.masked', 'BSgenome.Mmusculus.UCSC.mm10.masked', 'EnsDb.Hsapiens.v75', 'EnsDb.Hsapiens.v86', 'EnsDb.Mmusculus.v79'))" -e "devtools::install_github('databio/simpleCache')" -e "devtools::install_github('databio/GenomicDistributions')"

library(GenomicDistributions)

# move conf files
COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY shiny-server.sh /usr/bin/shiny-server.sh

# get the code
WORKDIR /srv/shiny-server
RUN git clone -b dev https://github.com/databio/LOLAweb.git

# add dir for cache
RUN mkdir LOLAweb/cache
RUN chown -R shiny:shiny LOLAweb/cache

# run the server setup script
CMD ["/usr/bin/shiny-server.sh"]
