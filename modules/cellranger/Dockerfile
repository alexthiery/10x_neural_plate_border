# Set the base image to Ubuntu 16.04
FROM ubuntu:16.04

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing cellranger"

# update apt-get and install wget
RUN  apt-get update \
  && apt-get install -y wget \
  && rm -rf /var/lib/apt/lists/*

# Install cellranger
RUN cd /opt/ && \
	wget -O cellranger-3.0.2.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-3.0.2.tar.gz?Expires=1595472309&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci0zLjAuMi50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE1OTU0NzIzMDl9fX1dfQ__&Signature=I6OgTvvSKl9mL0PoUeAoFuec9Ng6EMfu7JB03uZrna2q7V4vsHdatFNnrepyCKHF7pz6rdNbspSl1jL4dgkDy3y1K3yNTHWX5DQTQCUZNyQuEMai1QQIwB1B9o0q5Sz2LZDrqnFzdZfdVr9~6FS5IYvi4bKLHf4DTIyiawwiZri1eb93rjDcE6-j6LPwQFVafXr7xXdV5C7OlwIHYLAnPflJJ4WqzxtQU53gtdPqdH4oEYa3gWu5ax7UBx3mTlyJOtwrZZNOMvsUYAiww~m45gCirEWzzJ6aADGuWDGxXtkTMuc9v3Nh~Qfk-T1qmu6mq8jEpwvoJbWop66-IaSzTQ__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA" && \	
	tar -xzvf cellranger-3.0.2.tar.gz && \
	rm -f cellranger-3.0.2.tar.gz

# Set path
ENV PATH /opt/cellranger-3.0.2:$PATH

