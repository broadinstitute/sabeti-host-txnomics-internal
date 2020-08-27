FROM ubuntu:18.04

RUN apt-get update && apt-get -y install wget curl

RUN mkdir -p /opt/stringtie2 && cd /opt/stringtie2/ && \
    wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.1.4.Linux_x86_64.tar.gz && \
    tar xzf stringtie-2.1.4.Linux_x86_64.tar.gz

ENV PATH=${PATH}:/opt/stringtie2/stringtie-2.1.4.Linux_x86_64/