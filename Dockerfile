FROM rocker/shiny:latest
MAINTAINER Markus List <mlist@health.sdu.dk>

#install system packages
RUN apt-get update
RUN apt-get install -y libxml2-dev redis-server libcurl4-gnutls-dev libssl-dev 

#install R packages
ADD install.R install.R
RUN R -e "source('install.R')"

#copy shiny app to work-dir
WORKDIR /srv/
RUN mkdir hitseekr
ADD . hitseekr

#update shiny server conf and configure it to run hitseekr in single app mode
RUN sed -i 's/site_dir \/srv\/shiny-server;/app_dir \/srv\/hitseekr;/g' /etc/shiny-server/shiny-server.conf

#download additional database files if needed
#RUN wget https://www.dropbox.com/s/rr432s5tbc1gsby/compound.mapping.sqlite3?dl=0
#RUN wget https://www.dropbox.com/s/e19hjwqdpt3ajva/rnahybrid.sqlite3?dl=0
#RUN wget https://www.dropbox.com/s/2ip1awimawovt4r/stitch_hsa_protein_chemical_links_v4.0.sqlite3?dl=0