# Set the base image to Ubuntu 16.04
FROM rocker/tidyverse:3.6.3

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing Rstudio and getopt"

# update apt-get and install
RUN R -e "install.packages('getopt')"