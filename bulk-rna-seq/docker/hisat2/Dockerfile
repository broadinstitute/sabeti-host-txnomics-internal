FROM ubuntu:18.04

RUN apt-get update && apt-get install -y wget build-essential zlib1g-dev zip

RUN mkdir -p /opt/tools
WORKDIR /opt/tools

RUN wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
RUN tar xjf samtools-1.10.tar.bz2
WORKDIR /opt/tools/samtools-1.10
RUN ./configure --without-curses --disable-bz2 --disable-lzma
RUN make -j 2
RUN mv samtools ..

WORKDIR /opt/tools/
ENV PATH=/opt/tools:${PATH}
RUN rm -rf samtools-1.10  samtools-1.10.tar.bz2

RUN wget -O hisat2-2.2.0-Linux_x86_64.zip https://cloud.biohpc.swmed.edu/index.php/s/hisat2-220-Linux_x86_64/download
RUN unzip hisat2-2.2.0-Linux_x86_64.zip
ENV PATH=/opt/tools/hisat2-2.2.0/:${PATH}

