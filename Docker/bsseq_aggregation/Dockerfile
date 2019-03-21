FROM debian:jessie

MAINTAINER Divy Kangeyan <divyswar01@g.harvard.edu>

RUN apt-get update &&  apt-get install -y build-essential zlib1g-dev unzip libxml2-dev libssl-dev libcurl4-openssl-dev libmariadb-client-lgpl-dev pandoc

# Install R-3.5.0 (See https://cran.r-project.org/bin/linux/debian)
RUN apt-get install -y libssl-dev libcurl4-openssl-dev && \
apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF' && \
echo "deb http://cloud.r-project.org/bin/linux/debian jessie-cran35/" >> /etc/apt/sources.list && \
apt-get update && \
apt-get install -y --force-yes r-base




#Installing necessary R packages
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite(c('BSgenome', 'Biostrings','GenomeInfoDb'), ask=F)"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite(c('HDF5Array','bsseq'), ask=F)"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite(c('foreach','optparse'), ask=F)"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite(c('readr','stringr', 'devtools'), ask=F)"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite(c('annotatr'), ask=F)"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite(c('TxDb.Mmusculus.UCSC.mm10.knownGene'), ask=F)"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite(c('org.Mm.eg.db'), ask=F)"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite(c('TxDb.Hsapiens.UCSC.hg38.knownGene'), ask=F)"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite(c('TxDb.Hsapiens.UCSC.hg19.knownGene'), ask=F)"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite(c('org.Hs.eg.db'), ask=F)"
RUN Rscript -e "devtools::install_github('aryeelab/scmeth', ref = '9713607')" 
#RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite(c('scmeth'), ask=F)"

ADD ./bismark_to_bsseq_file_list.R bismark_to_bsseq_file_list.R

CMD ["/bin/bash"]