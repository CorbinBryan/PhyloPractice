FROM ubuntu

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \\
    & apt-get install mafft iqtree perl 

ENV USER root 

WORKDIR /root/ 

COPY AmanitaAccessionUSETHIS.txt /root 

RUN 
