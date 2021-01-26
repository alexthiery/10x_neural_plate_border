FROM rocker/tidyverse:3.6.3

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements to run 10x downstream analysis in Seurat and Antler"

# Install apt packages
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
 git=1:2.20.1-2+deb10u3 \
 apt-utils=1.8.2.1


 RUN R -e 