FROM  alexthiery/10x-npb-schelper:latest

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements to run 10x downstream analysis in Seurat"


RUN alias python=python3.8
RUN apt install python3-pip -y
RUN pip3 install --user phate
# Install phateR
RUN R -e "install.packages(c('phateR', 'Rmagic'))"
