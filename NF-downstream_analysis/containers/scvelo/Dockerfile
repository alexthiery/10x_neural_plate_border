FROM continuumio/miniconda3:4.6.14

ENV MPLCONFIGDIR=/tmp/matplotlib
ENV NUMBA_CACHE_DIR=/tmp/numbacache

# Install apt packages
RUN apt-get update \
      && apt-get install -y --no-install-recommends \
      apt-utils \
      build-essential \
      procps

# Create the environment:
COPY environment.yml .

RUN conda update -n base -c defaults conda && \
    conda env create -f environment.yml && \
    echo "conda activate scvelo" >> ~/.bashrc

RUN pip install git+https://github.com/theislab/cellrank@dev
    