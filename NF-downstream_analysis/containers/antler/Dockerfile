FROM alexthiery/10x-npb-seurat:base-1.5

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements to run 10x downstream analysis in Seurat and Antler"

# RUN   R -e "devtools::install_github('juliendelile/Antler', dependencies = TRUE)"

RUN   R -e "devtools::install_version('ArrayBin', version = '0.2')"
RUN   R -e "devtools::install_github('juliendelile/Antler', dependencies = TRUE)"
RUN   R -e "devtools::install_github('alexthiery/scHelper', dependencies = TRUE)"