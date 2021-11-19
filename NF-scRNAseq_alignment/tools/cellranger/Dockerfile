# Set the base image to Ubuntu 16.04
FROM ubuntu:16.04

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing cellranger"

# update apt-get and install
RUN  apt-get update \
  && apt-get install -y wget \
  && rm -rf /var/lib/apt/lists/*

# Install cellranger
RUN cd /opt/ && \
  wget -O cellranger-4.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-4.0.0.tar.gz?Expires=1618014446&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci00LjAuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2MTgwMTQ0NDZ9fX1dfQ__&Signature=To8mvt53qmeKTZA5gU2DMrukXjNXq8tAP1q2~Tv992hKerFjDJWYOuBc0vyWOT3IyI6Y36EnJHrNn5ni5nBDCOJad~MO0hbr1WZTDSApa77ml40xprY~J7hPDmfx0XX9Y8wgPwWcbfXdwFa-Rcjj2M~x9dpgQ4-aiPiW-OzOfyIhdrdvifp6g9IfUV527cMgeEYC1e3~neFQtWrfaO7FZLvtlpqSieQ04TjPb9gx9dDY2EpsI52fxhvKFnhF2j8Iev0qjfDGZmxhdxhIOaDPfG5w4yXlrbMJkEmGnIEA7XiL9kWxkiyUU7rU7st7dgIbOccpQxbPZIXhiYhrjVas2g__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA" && \	
	tar -xzvf cellranger-4.0.0.tar.gz && \
	rm -f cellranger-4.0.0.tar.gz

# Set path
ENV PATH /opt/cellranger-4.0.0:$PATH