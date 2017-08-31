
mkdir bismark_test
cd bismark_test

GENOME_FOLDER=/data/aryee/pub/genomes/bismark_index/grch38

SAMPLE_NAME=small_01
FASTQ_R1=../testdata/small_01_R1.fastq.gz
FASTQ_R2=../testdata/small_01_R2.fastq.gz

module load aryee/samtools-1.3.1
time /apps/lab/aryee/Bismark-0.18.2/bismark --genome $GENOME_FOLDER --basename $SAMPLE_NAME -1 $FASTQ_R1 -2 $FASTQ_R2

time /apps/lab/aryee/Bismark-0.18.2/bismark_methylation_extractor --gzip --bedGraph --buffer_size 4G --cytosine_report --genome_folder $GENOME_FOLDER ${SAMPLE_NAME}_pe.bam



library("BSgenome.Hsapiens.UCSC.hg38")
library("readr")
library("stringr")
library("HDF5Array")
library("bsseq")

dir="."
sample="tmp"


# Get phenodata
extract <- function(pattern, lines) {
	idx <- grep(pattern, lines)[1]
	line <- lines[idx]
	m <- str_match(line, pattern) # The extracted match is in column 2
	if (ncol(m)==2) {
		return(m[1,2])
	} else {
		return("")
	}
}

patterns <- c(  total_reads = "Sequence pairs analysed in total:\t(.*)",
				uniquely_aligned_reads = "Number of paired-end alignments with a unique best hit:\t(.*)",
				non_uniquely_aligned_reads = "Sequence pairs did not map uniquely:\t(.*)",
				OT_reads = "CT/GA/CT:\t(.*)\t",
				OB_reads = "CT/GA/GA:\t(.*)\t",
				CTOT_reads = "GA/CT/CT:\t(.*)\t",
				CTOB_reads = "GA/CT/GA:\t(.*)\t",
				aligner_version = "\\(version: (.*)\\)",
				aligner_options = "specified options: (.*)",
				genome = "bisulfite genome of (.*) with"
				)

bismark_report <- file.path(dir, paste0(sample, "_PE_report.txt"))
report_lines <- readLines(bismark_report)				
pd <- data.frame(row.names=sample,
			t(as.data.frame(sapply(patterns, extract, report_lines), stringsAsFactors=FALSE)))


# Get methylation matrices
cpg <- DNAString("CG")
system.time(cpg_gr <- vmatchPattern(cpg, Hsapiens))
cpg_gr <- keepStandardChromosomes(cpg_gr, pruning.mode= "coarse")
seqlevelsStyle(cpg_gr) <- "Ensembl"

# On the plus strand we keep the left-most position of the match
# On the minus strand we keep the right-most position of the match

s <- start(cpg_gr)
e <- end(cpg_gr)
plus_idx <- as.logical(strand(cpg_gr)=="+")
minus_idx <- as.logical(strand(cpg_gr)=="-")
e[plus_idx] <- s[plus_idx]  # Plus strand
s[minus_idx] <- e[minus_idx] # Minus strand
start(cpg_gr) <- s
end(cpg_gr) <- e

tab <- read_tsv(file.path(dir, paste0(sample, "_pe.bismark.cov.gz")), 
			col_types = "ciiiii", 
			col_names=c("chr", "pos", "pos2", "meth_percent", "m_count", "u_count"))
tab_gr <- GRanges(tab$chr, IRanges(tab$pos, tab$pos))
m <- rep(NA, length(cpg_gr))
cov <- rep(0, length(cpg_gr))
ovl <- findOverlaps(tab_gr, cpg_gr)
m[subjectHits(ovl)]   <- tab$m_count[queryHits(ovl)]
cov[subjectHits(ovl)] <- tab$m_count[queryHits(ovl)] + tab$u_count[queryHits(ovl)]





# Save bsseq object to rds using HDF5Arrays for on-disk storage
rda <- paste0(sample, ".rda")
hdf5_file <- paste0(sample, ".hdf5")
if(file.exists(rda)) file.remove(rda)
if(file.exists(hdf5_file)) file.remove(hdf5_file)

hdf5_m <- writeHDF5Array(matrix(m), file=hdf5_file)
hdf5_cov <- writeHDF5Array(matrix(cov), file=hdf5_file)
bs <- BSseq(gr=cpg_gr, M=hdf5_m, Cov=hdf5_cov, pData=pd, sampleNames=sample)
message("Saving bs to ", rda)
save(bs, file=rda)

message("Done.")






    # Get the phenotypic data
    readInfo<-read.table(readmetric_file,sep=' ')
    bsconvInfo<-read.table(bsConv_file,sep=' ', header=TRUE)
    pd<-data.frame(row.names=sample,totalReads=as.numeric(readInfo[2]),
                    mappedReads=as.numeric(readInfo[3]),
                    bsconversion=as.numeric(bsconvInfo[1,1]),
                    stringsAsFactors=FALSE)




time /apps/lab/aryee/Bismark-0.18.2/bismark_methylation_extractor --bedGraph --buffer_size 4G --genome_folder $GENOME_FOLDER tmp_pe.bam


time /apps/lab/aryee/Bismark-0.18.2/bismark --genome $GENOME_FOLDER -1 ../testdata/small_02_R1.fastq.gz -2 ../testdata/small_02_R2.fastq.gz
time /apps/lab/aryee/Bismark-0.18.2/bismark --genome $GENOME_FOLDER -1 ../testdata/small_03_R1.fastq.gz -2 ../testdata/small_03_R2.fastq.gz


bismark_methylation_extractor
