# Docker image based on Ubuntu
FROM ubuntu:latest
MAINTAINER Nate Matteson <natem@scripps.edu>

RUN apt-get update
RUN apt-get install -y python3-setuptools python3-docutils python3-flask default-jre gzip wget git python3-pip
RUN apt-get install -y bwa tabix

# Install samtools, hopefully
ENV SAMTOOLS_INSTALL_DIR=/opt/samtools
WORKDIR /tmp
RUN wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && \
  tar --bzip2 -xf samtools-1.3.1.tar.bz2
WORKDIR /tmp/samtools-1.3.1
RUN ./configure --enable-plugins --prefix=$SAMTOOLS_INSTALL_DIR && \
  make all all-htslib && \
  make install install-htslib
WORKDIR /
RUN ln -s $SAMTOOLS_INSTALL_DIR/bin/samtools /usr/bin/samtools && \
  rm -rf /tmp/samtools-1.3.1

ENV JAVA_HOME  /usr/lib/jvm/java-8-openjdk-amd64

# Clone the iVar repository.
WORKDIR /home/user
RUN mkdir /wd
RUN git clone https://github.com/watronfire/iVar.git
WORKDIR /home/user/iVar

# Install python requirements.
RUN pip3 install -r requirements.txt
RUN pip3 install snakemake

# Initiate docker with the command:
# docker run -it --name:iVar -v <pathToWD>:/home/user/wd watronfire/ivar