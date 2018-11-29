FROM debian:jessie

RUN apt-get install python

# Install procps (to get free and top)
RUN apt-get install -y procps

FROM python:3

RUN pip install pandas
RUN pip install numpy

ADD ./pdr.py pdr.py