FROM  rocker/tidyverse:4.0.3

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements to run 10x downstream analysis in Seurat and Antler"

# Install apt packages
RUN apt-get update \
      && apt-get install -y --no-install-recommends \
      apt-utils \
      libhdf5-dev=1.10.4+repack-10 \
      libboost-filesystem1.67-dev \
      libxtst6

RUN   R -e "install.packages('Seurat', dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
      R -e "install.packages('future', dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
      R -e "install.packages('cowplot', dependencies=TRUE, repos='http://cran.rstudio.com/')"  && \
      R -e "install.packages('clustree', dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
      R -e "install.packages('gridExtra', dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
      R -e "install.packages('getopt', dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
      R -e "install.packages('pheatmap', dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
      R -e "remotes::install_github('carmonalab/STACAS')" \
      R -e "BiocManager::install('limma')"