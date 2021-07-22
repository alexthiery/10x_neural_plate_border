FROM  alexthiery/10x-npb-seurat:latest

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements to run 10x downstream analysis in Seurat"

RUN   apt-get install python3-pip

RUN   pip install phate && \
      R -e "install.packages(c('phateR', 'Rmagic'))" && \
      R -e "devtools::install_github('alexthiery/scHelper', dependencies = TRUE)"