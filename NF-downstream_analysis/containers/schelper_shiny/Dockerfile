FROM  alexthiery/10x-npb-schelper:latest

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements to run 10x downstream analysis in Seurat and scHelper"

# Install scHelper
RUN   R -e "install.packages('shiny')" && \
      R -e "install.packages('shinydashboard')" && \
      R -e "install.packages('bs4Dash')"
