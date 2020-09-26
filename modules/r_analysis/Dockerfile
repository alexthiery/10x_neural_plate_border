FROM rocker/tidyverse:3.6.3

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements to run 10x downstream analysis in Seurat and Antler"

# Install apt packages
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
 git=1:2.20.1-2+deb10u3 \
 apt-utils=1.8.2.1 \
 unzip=6.0-23+deb10u1 \
 procps=2:3.3.15-2 \
 build-essential=12.6 \
 zlib1g-dev=1:1.2.11.dfsg-1 \
 libxt-dev \
 libhdf5-dev=1.10.4+repack-10 \
 libboost-filesystem1.67-dev \
 libmagick++-dev

RUN R -e "BiocManager::install(version = '3.10')"
RUN R -e "BiocManager::install('pcaMethods')"
RUN R -e "BiocManager::install('limma')"

RUN R -e "devtools::install_github('velocyto-team/velocyto.R')"
RUN R -e "devtools::install_github('juliendelile/Antler', dependencies = TRUE)"

RUN R -e "install.packages('Seurat')"
RUN R -e "install.packages('future')"
RUN R -e "install.packages('cowplot')"
RUN R -e "install.packages('clustree')"
RUN R -e "install.packages('gridExtra')"
RUN R -e "install.packages('getopt')"
RUN R -e "install.packages('pheatmap')"
RUN R -e "remotes::install_github('carmonalab/STACAS')"