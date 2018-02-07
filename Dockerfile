# Docker image based on Ubuntu
FROM ubuntu:16.04
MAINTAINER Nate Matteson <natem@scripps.edu>

RUN apt-get -qq update
RUN apt-get install -qqy python3-setuptools python3-docutils python3-flask build-essential zlib1g-dev libncurses5-dev default-jre gunzip git

ENV JAVA_HOME  /usr/lib/jvm/java-8-openjdk-amd64

# Clone the iVar repository.
WORKDIR /home
RUN git clone https://github.com/watronfire/iVar.git
WORKDIR /home/iVar

# Install python requirements.
RUN pip3 install -r requirements.txt