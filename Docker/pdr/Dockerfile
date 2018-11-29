FROM ubuntu

# Update
RUN apt-get update
RUN apt-get install -y python 
RUN apt-get install -y python-pip
 
# Install app dependencies
RUN pip install --upgrade pip

RUN pip install pandas
RUN pip install numpy

ADD ./pdr.py pdr.py