FROM continuumio/miniconda3:4.6.14

ENV MPLCONFIGDIR=/tmp/matplotlib
ENV NUMBA_CACHE_DIR=/tmp/numbacache

# Install packages and make conda env
RUN   apt-get update \
      && apt-get install -y --no-install-recommends \
      apt-utils \
      build-essential \
      procps && \
      rm -rf /var/lib/apt/lists/* && \
      conda create -n scvelo && \
      conda install --yes --freeze-installed -n scvelo python=3.8.5 scvelo numba=0.52.0 notebook ipywidgets --channel bioconda --channel defaults --channel conda-forge && \
      conda clean -afy && \
      echo "conda activate scvelo" >> ~/.bashrc