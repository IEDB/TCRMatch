FROM ubuntu:20.04

ENV PACKAGES python3 python3-pip git cmake

#avoid timezone graphical interaction
RUN ln -snf /usr/share/zoneinfo/$CONTAINER_TIMEZONE /etc/localtime && echo $CONTAINER_TIMEZONE > /etc/timezone

RUN apt-get update && \
    apt-get install -y ${PACKAGES} && \
    apt-get clean

RUN git clone https://github.com/IEDB/TCRMatch.git

RUN cd ./TCRMatch && cmake . && cmake --build .

