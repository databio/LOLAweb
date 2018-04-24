FROM databio/shinybase:latest

MAINTAINER VP Nagraj "vpnagraj@virginia.edu"

# move conf files
COPY shiny-server.conf  /etc/shiny-server/shiny-server.conf
COPY shiny-server.sh /usr/bin/shiny-server.sh

# get the code
WORKDIR /srv/shiny-server
RUN git clone https://github.com/databio/LOLAweb.git

# add dir for cache
RUN mkdir LOLAweb/cache
RUN chown -R shiny:shiny LOLAweb/cache

# run the server setup script
CMD ["/usr/bin/shiny-server.sh"]
