FROM python:3.7.4

RUN mkdir -p /opt/umi_collapser

WORKDIR /opt/umi_collapser

RUN git clone --branch 0.1.0 https://github.com/broadinstitute/umi-collapser-internal.git
RUN cd umi-collapser-internal && pip install .

ENV PATH=$PATH:/opt/umi_collapser/umi-collapser-internal/utils/