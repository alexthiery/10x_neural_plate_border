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
	wget -O cellranger-4.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-4.0.0.tar.gz?Expires=1599603640&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci00LjAuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE1OTk2MDM2NDB9fX1dfQ__&Signature=oNxXDS7Tlm5gqhaLsZ1ASfInlL-FoS3wmUHTKqAScdhgv54BrkYtEW1oVBGDAUiqrmOehgnYwRrGCk6E2~jeguxP5x1dmVYSrTZkxiTG1R7lfvOlivvPnTGn0AkqiXatzQOHngGPthboJuc-g35xkw-UUXetT5RDY9T6S4bHFnsR7TPgr4iMLAWC2WJLatTB3QN3OXAVjk2~WGUtMZihM7akWCXynqj3Gpkbh44R1Yz~7m2lvkI35qfdt~rUlJOmHysati72UdQsXMPITWWaBDeqcv0sHUItyEF~3YRnacrQQHrtBaMLgDFiRCIf1zHLnWrsqvYcaeqOm511plfr~A__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA" && \	
	tar -xzvf cellranger-4.0.0.tar.gz && \
	rm -f cellranger-4.0.0.tar.gz

# Set path
ENV PATH /opt/cellranger-4.0.0:$PATH