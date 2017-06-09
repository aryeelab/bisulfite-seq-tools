FROM debian:jessie
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
RUN mkdir /Rscripts/
COPY *.R /Rscripts/

#Installing necessary R packages
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('bsseq');"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('BiocGenerics');"
RUN Rscript -e "install.packages('optparse', repos='http://R-Forge.R-project.org')"
RUN Rscript -e "install.packages('ggplot2',repos='http://R-Forge.R-project.org')"
RUN Rscript -e "install.packages('foreach',repos='http://R-Forge.R-project.org')"
RUN Rscript -e "install.packages('yaml',repos='https://mirror.las.iastate.edu/CRAN/')"
RUN Rscript -e "install.packages('devtools',repos='https://mirror.las.iastate.edu/CRAN/')"
RUN apt-get install -y libxml2-dev
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('XML')"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('genefilter')"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('rtracklayer')"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('BSgenome')"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('Ringo')"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('Repitools')"
RUN mkdir /executable_files/
COPY collect_bsconv_metrics.sh /executable_files/
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('knitr')"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('AnnotationHub')"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('BSgenome.Mmusculus.UCSC.mm10')"
RUN Rscript -e "install.packages('rmarkdown',repos='https://mirror.las.iastate.edu/CRAN/')"
COPY scmeth-master /Rscripts/scmeth
RUN  R CMD INSTALL /Rscripts/scmeth
COPY Rscripts/*.R /Rscripts/
RUN apt-get install -y pandoc
