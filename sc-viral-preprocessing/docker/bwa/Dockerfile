FROM ubuntu:18.04

RUN apt-get update && apt-get install -y default-jre wget samtools build-essential zlib1g-dev

RUN mkdir -p /opt/bwa

WORKDIR /opt/bwa

RUN wget "https://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fbio-bwa%2Ffiles%2Fbwa-0.7.17.tar.bz2%2Fdownload%3Fuse_mirror%3Dcfhcable%26r%3Dhttps%253A%252F%252Fsourceforge.net%252Fprojects%252Fbio-bwa%252Ffiles%252Fbwa-0.7.17.tar.bz2%252Fdownload&ts=1595519387" -O bwa-0.7.17.tar.bz2

RUN tar xvjf bwa-0.7.17.tar.bz2 && cd bwa-0.7.17 && make && mv bwa ..

ENV PATH=$PATH:/opt/bwa