FROM ubuntu:18.04

RUN apt-get update && apt-get -y install wget
RUN wget "https://sourceforge.net/projects/subread/files/subread-2.0.1/subread-2.0.1-Linux-x86_64.tar.gz/download" -O "subread-2.0.1-Linux-x86_64.tar.gz"
RUN tar xvzf subread-2.0.1-Linux-x86_64.tar.gz
RUN mkdir -p /opt/subread/
RUN mv subread-2.0.1-Linux-x86_64/bin/* /opt/subread/
ENV PATH=${PATH}:/opt/subread/