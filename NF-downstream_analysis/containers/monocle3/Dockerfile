FROM  alexthiery/10x-npb-schelper:latest

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements to run 10x downstream analysis in Seurat"

RUN   apt-get update \
      && apt-get install -y --no-install-recommends \
      libboost1.71-dev \
      libboost-dev \
	libkml-dev \
      libudunits2-dev \
      libgdal-dev \
      libgeos-dev \
      libproj-dev \
      libharfbuzz-dev \
      libfribidi-dev


RUN   conda install -c anaconda pip && \
      pip install phate && \
      R -e "BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats', 'limma', 'S4Vectors', 'SingleCellExperiment', \
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))" && \
      R -e "devtools::install_github('cole-trapnell-lab/leidenbase')" && \
      R -e "devtools::install_github('cole-trapnell-lab/monocle3')"
      