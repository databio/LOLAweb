FROM databio/shinybase:latest

MAINTAINER VP Nagraj "nagraj@nagraj.net"

COPY run.sh /usr/bin/run.sh

# get the code
WORKDIR /srv/shiny-server
RUN git clone --branch shinyqueue https://github.com/databio/LOLAweb.git

# add dir for tempdir
RUN mkdir LOLAweb/temp
RUN chown -R shiny:shiny LOLAweb/temp

CMD ["/usr/bin/run.sh"]
