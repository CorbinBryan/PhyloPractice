FROM openjdk:8

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y mafft wget perl 

ENV USER root 

WORKDIR /root/ 

COPY ./wi_amanita.fasta /root 

# bro im lazy, precombiled binaries it is
RUN wget https://github.com/iqtree/iqtree2/releases/download/v2.2.2.6/iqtree-2.2.2.6-Linux.tar.gz \ 
    && mkdir ./iqtree \
    && tar -xzf ./iqtree-2.2.2.6-Linux.tar.gz 

ENV PATH=/root/iqtree-2.2.2.6-Linux/bin/:$PATH