FROM debian:jessie
MAINTAINER Andrew Dunford <adunford@broadinstitute.org>


RUN apt-get update &&  apt-get install  -y build-essential zlib1g-dev unzip libxml2-dev

# Install R-3.4.0 (See https://cran.r-project.org/bin/linux/debian)
RUN apt-get install -y libssl-dev libcurl4-openssl-dev && \
	apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF' && \
	echo "deb http://cloud.r-project.org/bin/linux/debian jessie-cran34/" >> /etc/apt/sources.list && \
	apt-get update && \
	apt-get install -y --force-yes r-base

#Installing necessary R packages

RUN Rscript -e "source('http://bioconductor.org/biocLite.R');required_libraries <- c('BSgenome','readr','stringr','HDF5Array','bsseq','foreach','optparse','Biostrings','GenomeInfoDb','devtools');sapply(required_libraries, biocLite, ask=F)"
RUN     Rscript -e "devtools::install_github('aryeelab/scmeth', ref = 'bismark')" && \
      apt-get install -y pandoc
ADD ./bismark_to_bsseq_file_list.R bismark_to_bsseq_file_list.R
