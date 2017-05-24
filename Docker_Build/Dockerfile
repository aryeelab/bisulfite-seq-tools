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

#Annoying Library stuff needed to install R-3.3.0
RUN apt-get install -y libcurl4-openssl-dev wget gfortran
WORKDIR /src
RUN wget http://www.bzip.org/1.0.6/bzip2-1.0.6.tar.gz
RUN mkdir -p /usr/local/packages6/libs/gcc/4.4.7/bzip2/1.0.6
RUN tar -xvzf ./bzip2-1.0.6.tar.gz
WORKDIR bzip2-1.0.6
RUN make -f Makefile-libbz2_so
RUN make
RUN make install PREFIX=/usr/local/packages6/libs/gcc/4.4.7/bzip2/1.0.6
RUN mv *.so* /usr/local/packages6/libs/gcc/4.4.7/bzip2/1.0.6/lib/

COPY bashrc_addendum /root/.bashrc
#RUN cat bashrc_addendum >> /root/.bashrc

#installing R-3.3.0
RUN mkdir -p /src/R/
WORKDIR /src/R/
RUN  wget  https://cran.r-project.org/src/base/R-3/R-3.3.0.tar.gz
RUN tar xzvf R-3.3.0.tar.gz
WORKDIR R-3.3.0
#RUN source ~/.bashrc && echo $LD_LIBRARY_PATH
#RUN source ~/.bashrc && echo $PATH
#RUN ./configure  --with-readline=no --with-x=no --libdir=/usr/local/packages6/libs/gcc/4.4.7/bzip2/1.0.6/lib/ --enable-R-shlib --with-blas --with-lapack
COPY install_R.sh .
RUN sh install_R.sh