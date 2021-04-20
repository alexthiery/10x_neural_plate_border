FROM  alexthiery/10x-npb-schelper:latest

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements to run 10x downstream analysis in Seurat"

RUN   pip install phate && \
      R -e "install.packages(c('phateR', 'Rmagic'))"