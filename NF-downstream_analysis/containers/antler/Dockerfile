# FROM alexthiery/10x-npb-seurat:base-1.5

# LABEL authors="alex.thiery@crick.ac.uk" \
#       description="Docker image containing all requirements to run 10x downstream analysis in Seurat and Antler"

# RUN   R -e "devtools::install_github('juliendelile/Antler', dependencies = TRUE)"

FROM  rocker/tidyverse:3.6.3

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

RUN   R -e "devtools::install_github('juliendelile/Antler', dependencies = TRUE)"