FROM ubuntu:18.04

RUN apt-get update && apt-get install -y default-jre wget samtools

RUN mkdir -p /opt/picard

WORKDIR /opt/picard

RUN wget https://github.com/broadinstitute/picard/releases/download/2.23.3/picard.jar -O picard.jar

ENV PICARD_JAR_PATH=/opt/picard/picard.jar