FROM  rocker/tidyverse:4.0.3

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements to run 10x downstream analysis in Seurat and Antler"

# Install apt packages
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
 git=1:2.20.1-2+deb10u3 \
 apt-utils=1.8.2.2

RUN   R -e "install.packages('Seurat', dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
      R -e "install.packages('future', dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
      R -e "install.packages('cowplot', dependencies=TRUE, repos='http://cran.rstudio.com/')"  && \
      R -e "install.packages('clustree', dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
      R -e "install.packages('gridExtra', dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
      R -e "install.packages('getopt', dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
      R -e "install.packages('pheatmap', dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
      R -e "remotes::install_github('carmonalab/STACAS')"