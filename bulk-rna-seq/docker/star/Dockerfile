FROM ubuntu:18.04

RUN apt-get update && apt-get -y install wget libgomp1
WORKDIR /opt/
RUN wget https://github.com/alexdobin/STAR/archive/2.7.5a.tar.gz && tar xvzf 2.7.5a.tar.gz && rm 2.7.5a.tar.gz
ENV PATH=/opt/STAR-2.7.5a/bin/Linux_x86_64/:${PATH}