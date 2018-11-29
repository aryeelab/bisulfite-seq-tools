FROM debian:jessie

RUN apt-get update &&  apt-get install  -y python


RUN apt-get install python
# Install procps (to get free and top)
RUN apt-get install -y procps

RUN pip install pandas
RUN pip install numpy

ADD ./pdr.py pdr.py