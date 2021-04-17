FROM  alexthiery/10x-npb-schelper:latest

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements to run 10x downstream analysis in Seurat"

RUN   /root/miniconda3/bin/python3.8 -m pip install phate && \
      R -e "install.packages('phateR')"