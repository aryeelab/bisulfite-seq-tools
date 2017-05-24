FROM debian
MAINTAINER Andrew Dunford <adunford@broadinstitute.org>


#install bsmap
RUN apt-get update &&  apt-get install  -y build-essential zlib1g-dev samtools git
RUN mkdir /src
WORKDIR /src
COPY bsmap-2.90.tgz .
RUN tar zvxf bsmap-2.90.tgz
WORKDIR bsmap-2.90

RUN make
RUN make install

#install MethyDackel
#COPY MethylDackel /src/MethylDackel
WORKDIR /src
RUN git clone https://github.com/dpryan79/MethylDackel.git && \
	cd MethylDackel && \
	git reset --hard ea29c16 && \
	make && \
	make install

# Install R-3.4.0 (See https://cran.r-project.org/bin/linux/debian)
RUN apt-get install -y libssl-dev libcurl4-openssl-dev
RUN apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF' && \
	echo "deb http://cloud.r-project.org/bin/linux/debian jessie-cran34/" >> /etc/apt/sources.list && \
	apt-get update && \
	apt-get install -y --force-yes r-base=3.4.0-1~jessiecran.0
