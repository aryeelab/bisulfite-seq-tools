FROM debian:jessie

# Install Java
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive \
    apt-get -y install default-jre-headless && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
    
RUN apt-get update &&  apt-get install  -y build-essential zlib1g-dev libncurses5-dev libncursesw5-dev git wget unzip python-pip python-dev

# Install cromwell    
RUN wget https://github.com/broadinstitute/cromwell/releases/download/30.2/cromwell-30.2.jar