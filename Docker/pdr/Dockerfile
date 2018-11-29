FROM debian:jessie

RUN apt-get update &&  apt-get install  -y python


RUN apt-get install python
# Install procps (to get free and top)
RUN apt-get install -y procps


ADD ./pdr.py pdr.py