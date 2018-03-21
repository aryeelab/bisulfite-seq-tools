FROM rocker/r-ver:devel

MAINTAINER Divy Kangeyan <divyswar01@g.harvard.edu>

RUN apt-get update &&  apt-get install -y build-essential zlib1g-dev unzip libxml2-dev libssl-dev libcurl4-openssl-dev libmariadb-client-lgpl-dev pandoc

#Installing necessary R packages
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite(c('BSgenome', 'Biostrings','GenomeInfoDb'), ask=F)"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite(c('HDF5Array','bsseq'), ask=F)"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite(c('foreach','optparse'), ask=F)"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite(c('readr','stringr', 'devtools'), ask=F)"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite(c('annotatr'), ask=F)"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite(c('TxDb.Mmusculus.UCSC.mm10.knownGene'), ask=F)"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite(c('org.Mm.eg.db'), ask=F)"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite(c('TxDb.Hsapiens.UCSC.hg38.knownGene'), ask=F)"
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite(c('org.Hs.eg.db'), ask=F)"
RUN Rscript -e "devtools::install_github('aryeelab/scmeth', ref = '312b72b')" 

ADD ./bismark_to_bsseq_file_list.R bismark_to_bsseq_file_list.R

CMD ["/bin/bash"]