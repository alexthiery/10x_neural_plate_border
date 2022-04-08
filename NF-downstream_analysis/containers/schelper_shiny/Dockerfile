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

# Add conda to $PATH
ENV PATH /opt/conda/bin:$PATH

# Install miniconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-py38_4.9.2-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# Update conda and install phate
RUN   conda update -n base -c defaults conda

# Set default python for reticulate
RUN   echo 'RETICULATE_PYTHON = /opt/conda/bin/python' >> /usr/local/lib/R/etc/Renviron

# Set CRAN and BioC repositories for package installations
RUN   echo 'local({library(BiocManager); r <- repositories(); r["CRAN"] <- "https://cran.rstudio.com/"; options(repos = r)}); \
            options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version$platform, R.version$arch, R.version$os)))' >> /usr/local/lib/R/etc/Rprofile.site
            
# Install scHelper
RUN   R -e "install.packages('tidyverse')" && \
      R -e "install.packages('shiny')" && \
      R -e "install.packages('bs4Dash')" && \
      R -e "install.packages('Seurat')" && \
      R -e "install.packages('mgcv')" && \
      R -e "install.packages('viridis')" && \
      R -e "install.packages('devtools')" && \
      R -e "devtools::install_github('alexthiery/scHelper@v0.2.4', dependencies = TRUE)"