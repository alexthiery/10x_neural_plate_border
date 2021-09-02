FROM  alexthiery/10x-npb-seurat:latest

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements to run 10x downstream analysis in Seurat and scHelper"

# Install scHelper
RUN   R -e "devtools::install_github('alexthiery/scHelper@v0.1.0', dependencies = TRUE)"
