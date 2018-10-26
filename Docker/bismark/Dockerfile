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


# Install R-3.4.0 (See https://cran.r-project.org/bin/linux/debian)
RUN apt-get install -y libssl-dev libcurl4-openssl-dev && \
	apt-key adv --keyserver keys.gnupg.net --recv-key 'E19F5F87128899B192B1A2C2AD5F960A256A04AF' && \
	echo "deb http://cloud.r-project.org/bin/linux/debian jessie-cran34/" >> /etc/apt/sources.list && \
	apt-get update && \
	apt-get install -y --force-yes r-base

#Installing necessary R packages
RUN Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('bsseq');" && \
    Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('BiocGenerics');" && \
    Rscript -e "install.packages('optparse', repos='http://R-Forge.R-project.org')" && \
    Rscript -e "install.packages('ggplot2',repos='http://R-Forge.R-project.org')" && \
    Rscript -e "install.packages('foreach',repos='http://R-Forge.R-project.org')" && \
    Rscript -e "install.packages('yaml',repos='https://mirror.las.iastate.edu/CRAN/')" && \
    Rscript -e "install.packages('devtools',repos='https://mirror.las.iastate.edu/CRAN/')" && \
    apt-get install -y libxml2-dev && \
    Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('XML')" && \
    Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('genefilter')" && \
    Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('rtracklayer')" && \
    Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('BSgenome')" && \
    Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('Ringo')" && \
    Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('Repitools')" && \
    Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('knitr')" && \
    Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('AnnotationHub')" && \
    Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('BSgenome.Mmusculus.UCSC.mm10')" && \
    Rscript -e "install.packages('rmarkdown',repos='https://mirror.las.iastate.edu/CRAN/')" && \
    Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('annotatr')" && \
    Rscript -e "source('http://bioconductor.org/biocLite.R');biocLite('scmeth')" 


# Install bedtools
RUN apt-get install -y python && \
    cd /src && \
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
RUN pip install cutadapt

# Install FastQC
RUN apt-get -y install openjdk-7-jre && \
    cd /src && \
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip && \
    unzip fastqc_v0.11.7.zip && \
    ln -s /src/FastQC/fastqc /usr/local/bin/

# Install Trimgalore
RUN cd /tmp && \
    wget https://github.com/FelixKrueger/TrimGalore/archive/0.4.5.zip && \
    unzip 0.4.5.zip && \
    cp TrimGalore-0.4.5/trim_galore /usr/local/bin && \
    chmod +x /usr/local/bin/trim_galore

# added pooling report script
ADD ./pool_report_files.R pool_report_files.R