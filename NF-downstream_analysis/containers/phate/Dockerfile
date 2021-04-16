FROM  alexthiery/10x-npb-schelper:latest

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements to run 10x downstream analysis in Seurat"

# Initialise reticulate compatible conda
RUN   R -e "reticulate::install_miniconda()" && \
      R -e "reticulate::py_install('phate', pip=TRUE)" && \
      R -e "install.packages('phateR')"