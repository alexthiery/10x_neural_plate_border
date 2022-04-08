FROM  rocker/tidyverse:4.0.4

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements to run 10x downstream analysis in Seurat"

# Install apt packages
RUN apt-get update \
      && apt-get install -y --no-install-recommends \
      apt-utils \
      build-essential \
      libhdf5-dev \
      libboost-filesystem1.67-dev \
      libxtst6 \
      libxt-dev
      # libfontconfig1-dev \
      # libfribidi-dev \
      # libxml2

# Install scHelper
RUN   R -e "install.packages('tidyverse')" && \
      R -e "install.packages('shiny')" && \
      R -e "install.packages('bs4Dash')" && \
      R -e "install.packages('Seurat')" && \
      R -e "install.packages('mgcv')" && \
      R -e "install.packages('viridis')" && \
      R -e "install.packages('devtools')" && \
      R -e "devtools::install_github('alexthiery/scHelper@v0.2.4', dependencies = TRUE)"