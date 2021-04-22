FROM  alexthiery/10x-npb-seurat:latest

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing packages to run scVelo from Seurat object"

RUN pip install -U scvelo && \
    pip install python-igraph louvain && \
    pip install pybind11 hnswlib && \
    conda install notebook

RUN R -e "remotes::install_github('mojaveazure/seurat-disk')" && \
    R -e "install.packages('SeuratWrappers')" && \
    R -e "remotes::install_github('satijalab/seurat-wrappers')"