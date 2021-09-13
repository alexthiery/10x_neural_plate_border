FROM ubuntu

RUN apt-get update \
      && apt-get install -y --no-install-recommends \
      apt-utils \
      build-essential \
      python3 \
      pip

RUN echo "alias python=python3" >> ~/.bashrc

RUN   pip install cmake==3.18.4 && \
      pip install numba==0.51.2 && \
      pip install PhenoGraph && \
      pip install palantir && \
      pip install jupyterlab