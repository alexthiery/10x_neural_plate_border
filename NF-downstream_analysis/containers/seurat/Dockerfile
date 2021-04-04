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
      libxt-dev \
      libfontconfig1-dev \
      libfribidi-dev

# Set CRAN and BioC repositories for package installations
RUN   echo 'local({library(BiocManager); r <- repositories(); r["CRAN"] <- "https://cran.rstudio.com/"; options(repos = r)}); \
            options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version$platform, R.version$arch, R.version$os)))' >> /usr/local/lib/R/etc/Rprofile.site

RUN   R -e "install.packages('future', dependencies=TRUE)" && \
      R -e "install.packages('cowplot', dependencies=TRUE)"  && \
      R -e "install.packages('clustree', dependencies=TRUE)" && \
      R -e "install.packages('gridExtra', dependencies=TRUE)" && \
      R -e "install.packages('getopt', dependencies=TRUE)" && \
      R -e "install.packages('pheatmap', dependencies=TRUE)" && \
      R -e "install.packages('Seurat', dependencies=TRUE)" && \
      # install STACAS for time being
      R -e "remotes::install_github('carmonalab/STACAS')"

RUN   R -e "devtools::install_version('ArrayBin', version = '0.2')" && \
      R -e "devtools::install_github('juliendelile/Antler', dependencies = TRUE)"