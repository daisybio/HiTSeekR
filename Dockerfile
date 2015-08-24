FROM rocker/shiny:latest
MAINTAINER Markus List <mlist@health.sdu.dk>

#copy shiny app to work-dir
WORKDIR /srv/shiny-server/
RUN mkdir RNAice
COPY . /RNAice

#install system packages
RUN sudo apt-get install libxml2-dev redis-server 

#install R packages
RUN R -e "source('RNAice/install.R')"
