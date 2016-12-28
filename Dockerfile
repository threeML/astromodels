FROM ubuntu:16.04

MAINTAINER Giacomo Vianello <giacomov@stanford.edu>

# Override the default shell (sh) with bash
RUN rm /bin/sh && ln -s /bin/bash /bin/sh

# Update repositories and install needed packages
# I use one long line so that I can remove all .deb files at the end
# before Docker writes out this layer

RUN apt-get update && apt-get install -y python2.7 python2.7-dev python-dev wget python-pip git python-tk libreadline6-dev libncurses5-dev xorg-dev gcc g++ gfortran perl-modules && apt-get clean 

RUN mkdir -p /heasoft/src

COPY heasoft-6.19src.tar.gz /heasoft/src

RUN cd /heasoft && mkdir build && cd src && tar zxvf heasoft-6.19src.tar.gz && cd heasoft-6.19/BUILD_DIR/ && ./configure --prefix=/heasoft/build && make && make install && rm -rf /heasoft/src

RUN pip install virtualenv

RUN apt-get install sudo

WORKDIR /
