FROM ubuntu:18.04

ENV PACKAGES python3-dev python3-pip

RUN apt-get update && \
    apt-get install -y ${PACKAGES} && \
    apt-get clean

#RUN pip3 install tcrmatch
