FROM ubuntu:20.04

ENV TZ=America/Los_Angeles
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update
RUN apt-get install -y git cmake g++ curl

RUN git clone https://github.com/IEDB/TCRMatch.git
RUN cd /TCRMatch/ && make
RUN cd /TCRMatch/scripts && ./update.sh
