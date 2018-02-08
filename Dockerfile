# Docker image based on Ubuntu
FROM ubuntu:latest
MAINTAINER Nate Matteson <natem@scripps.edu>

RUN apt-get update
RUN apt-get install -y python3-setuptools python3-docutils python3-flask default-jre gzip git python3-pip

ENV JAVA_HOME  /usr/lib/jvm/java-8-openjdk-amd64

# Clone the iVar repository.
WORKDIR /home/user
RUN git clone https://github.com/watronfire/iVar.git
WORKDIR /home/user/iVar

# Install python requirements.
RUN pip3 install -r requirements.txt
RUN pip3 install snakemake