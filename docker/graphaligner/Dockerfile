FROM ubuntu:22.04

MAINTAINER rlorigro@broadinstitute.edu

ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

## Foundational
RUN apt-get update
RUN apt-get install -y --reinstall ca-certificates
RUN update-ca-certificates -f

## Basic software
RUN apt-get install -y --no-install-recommends wget
RUN apt-get install -y --no-install-recommends curl
RUN apt-get install -y --no-install-recommends cmake
RUN apt-get install -y --no-install-recommends autoconf
RUN apt-get install -y --no-install-recommends build-essential
RUN apt-get install -y --no-install-recommends bzip2
RUN apt-get install -y --no-install-recommends git
RUN apt-get install -y --no-install-recommends sudo
RUN apt-get install -y --no-install-recommends pkg-config
RUN apt-get install -y --no-install-recommends zlib1g-dev
RUN apt-get install -y --no-install-recommends libcurl4-openssl-dev
RUN apt-get install -y --no-install-recommends libbz2-dev
RUN apt-get install -y --no-install-recommends libncurses5-dev
RUN apt-get install -y --no-install-recommends liblzma-dev

## Add time library for basic profiling
RUN apt-get install -y --no-install-recommends time

## Install python3
RUN apt-get install -y --no-install-recommends python3
RUN apt-get install -y --no-install-recommends python3-pip

WORKDIR /software

## Install miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py310_23.3.1-0-Linux-x86_64.sh -O miniconda.sh
RUN bash miniconda.sh -b -f -p /software/miniconda
ENV PATH=$PATH:/software/miniconda/bin

## Install graphaligner with conda
RUN conda install -c bioconda graphaligner

## Force rebuild with time sensitive command
ADD http://date.jsontest.com /etc/builddate

## Copy repository python scripts into software dir
ADD ../../scripts /software

WORKDIR /software

ENTRYPOINT ["GraphAligner"]

