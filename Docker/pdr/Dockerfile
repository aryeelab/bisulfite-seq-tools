FROM debian:jessie

RUN apt-get update &&  apt-get install  -y build-essential zlib1g-dev libncurses5-dev libncursesw5-dev git wget unzip python-pip python-dev

# Install procps (to get free and top)
RUN apt-get install -y procps

# Install Bowtie2
RUN apt-get install -y libtbb-dev && \
	mkdir /src && \
	cd /src && \
	wget -O bowtie2.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.2/bowtie2-2.3.2-linux-x86_64.zip/download && \
	unzip bowtie2.zip && \
	ln -s /src/bowtie2-2.3.2/bowtie2* /usr/local/bin

# Install Bismark
RUN cd /src && \
	wget https://github.com/FelixKrueger/Bismark/archive/0.18.2.tar.gz && \
	tar zxf 0.18.2.tar.gz && \
	ln -s /src/Bismark-0.18.2/bismark* /usr/local/bin/

ADD ./pdr.py pdr.py