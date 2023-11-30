FROM openjdk:8

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    & apt-get install mafft iqtree wget perl 

ENV USER root 

WORKDIR /root/ 

COPY AmanitaAccessionUSETHIS.txt /root 

RUN yes | sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"

RUN wget https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/xtract.Linux.gz \
    && gzip -d xtract.Linux.gz \
    && chmod +x xtract.Linux \
    && mv xtract.Linux ~/./bin 

ENV PATH=/root/edirect:$PATH

# bro im lazy, precombiled binaries it is
RUN wget https://github.com/iqtree/iqtree2/releases/download/v2.2.2.6/iqtree-2.2.2.6-Linux.tar.gz \ 
    && mkdir ./iqtree \
    && tar -xzf ./iqtree-2.2.2.6-Linux.tar.gz 

ENV PATH=/root/iqtree-2.2.2.6-Linux/bin/:$PATH