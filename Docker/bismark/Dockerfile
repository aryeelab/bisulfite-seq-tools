FROM ubuntu:20.04

ARG VERSION=__specify_in_build_arg__

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \ 
    apt-get install -y build-essential zlib1g-dev \
                       libxml2-dev libncurses5-dev libncursesw5-dev git wget unzip \
                       python python-dev python3-pip \
                       libssl-dev libcurl4-openssl-dev dirmngr gnupg \
                       apt-transport-https ca-certificates software-properties-common\
                       fastqc

# Install procps (to get free and top)
RUN apt-get install -y procps

# Install Bowtie2
RUN apt-get install -y libtbb-dev && \
	mkdir /src && \
	cd /src && \
	wget -O bowtie2.zip https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/bowtie2-2.4.2-linux-x86_64.zip/download && \
	unzip bowtie2.zip && \
	ln -s /src/bowtie2-2.4.2/bowtie2* /usr/local/bin

# Install R
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
	add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' && \
	apt-get update && \
	apt-get install -y --force-yes r-base
	
#Installing necessary R packages
RUN Rscript -e "install.packages('BiocManager');" && \
    Rscript -e "install.packages('optparse', repos='http://R-Forge.R-project.org')" && \
    Rscript -e "install.packages('ggplot2',repos='http://R-Forge.R-project.org')" && \
    Rscript -e "install.packages('foreach',repos='http://R-Forge.R-project.org')" && \
    Rscript -e "install.packages('yaml',repos='https://mirror.las.iastate.edu/CRAN/')" && \
    Rscript -e "install.packages('devtools',repos='https://mirror.las.iastate.edu/CRAN/')" && \
    Rscript -e "install.packages('rmarkdown',repos='https://mirror.las.iastate.edu/CRAN/')" && \
    Rscript -e "BiocManager::install('bsseq');" && \
    Rscript -e "BiocManager::install('BiocGenerics');" && \
    Rscript -e "BiocManager::install('XML')" && \
    Rscript -e "BiocManager::install('genefilter')" && \
    Rscript -e "BiocManager::install('rtracklayer')" && \
    Rscript -e "BiocManager::install('BSgenome')" && \
    Rscript -e "BiocManager::install('Ringo')" && \
    Rscript -e "BiocManager::install('Repitools')" && \
    Rscript -e "BiocManager::install('knitr')" && \
    Rscript -e "BiocManager::install('AnnotationHub')" && \
    Rscript -e "BiocManager::install('BSgenome.Mmusculus.UCSC.mm10')" && \
    Rscript -e "BiocManager::install('annotatr')" && \
    Rscript -e "BiocManager::install('scmeth')" 

# Install bedtools
RUN cd /src && \
    wget https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz && \
    tar -zxf bedtools-2.26.0.tar.gz && \
    cd bedtools2 && \
    make && \
	ln -s /src/bedtools2/bin/bedtools /usr/local/bin

# Install samtools
RUN cd /src && \
    git clone https://github.com/samtools/htslib.git && \
    cd htslib && \
    git checkout tags/1.3 && \
    make && \
    cd .. && \
    git clone https://github.com/samtools/samtools.git && \
    cd samtools && \
    git checkout tags/1.3 && \
    make && \
    ln -s /src/samtools/samtools /usr/local/bin

# Install bedgraphToBigwig
RUN cd /usr/local/bin && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig && \
    chmod +x bedGraphToBigWig

# Install cutadapt (for Trimgalore)
RUN pip3 install cutadapt

# Install Trimgalore
RUN cd /tmp && \
    wget https://github.com/FelixKrueger/TrimGalore/archive/0.6.5.zip && \
    unzip 0.6.5.zip && \
    cp TrimGalore-0.6.5/trim_galore /usr/local/bin && \
    chmod +x /usr/local/bin/trim_galore

# Install Bismark
RUN cd /src && \
	wget https://github.com/FelixKrueger/Bismark/archive/0.22.3.tar.gz && \
	tar zxf 0.22.3.tar.gz && \
	ln -s /src/Bismark-0.22.3/bismark* /usr/local/bin/

# added pooling report script
ADD ./pool_report_files.R pool_report_files.R

RUN echo $VERSION > /VERSION
