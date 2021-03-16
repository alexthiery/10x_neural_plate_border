FROM  rocker/tidyverse:4.0.0

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements to run 10x downstream analysis in Seurat and Antler"

COPY ./.Rprofile /usr/local/lib/R/etc/Rprofile.site
# COPY ./.Rprofile /home/rstudio

# ENV R_PROFILE=/usr/local/lib/R/.Rprofile


# Install apt packages
RUN apt-get update \
      && apt-get install -y --no-install-recommends \
      apt-utils \
      build-essential \
      libhdf5-dev \
      libboost-filesystem1.67-dev \
      libxtst6 \
      libxt-dev \
      libfreetype6-dev




RUN   R -e "install.packages('Seurat', dependencies=TRUE)" && \
      R -e "install.packages('future', dependencies=TRUE)" && \
      R -e "install.packages('cowplot', dependencies=TRUE)"  && \
      R -e "install.packages('clustree', dependencies=TRUE)" && \
      R -e "install.packages('gridExtra', dependencies=TRUE)" && \
      R -e "install.packages('getopt', dependencies=TRUE)" && \
      R -e "install.packages('pheatmap', dependencies=TRUE)" && \
      R -e "remotes::install_github('carmonalab/STACAS')"