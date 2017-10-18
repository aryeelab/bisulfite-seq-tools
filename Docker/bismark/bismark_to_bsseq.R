# This script combines multiple bisulfite sequencing samples processed by the Bismark aligner
# into a single BSseq object saved as an HDF5SummarizedExperiment.
#
# Usage example
# Rscript --vanilla bismark_to_bsseq.R -i small_01,small_02,small_03 -o test_bs

library("BSgenome.Hsapiens.UCSC.hg38")
library("readr")
library("stringr")
library("HDF5Array")
library("bsseq")
library("foreach")
#getOption("DelayedArray.block.size")
options(DelayedArray.block.size=4500e6)
library("optparse")

option_list = list(
  make_option(c("-i","--input_samples"),type='character',default=NULL,
              help="Output from bismark_methylation_extractor. Use comma-separated sample names as passed to bismark with the --basename option"),
  make_option(c("-o","--output_bsseq_dir"),type='character',default=NULL,
              help="Directory in which to save the BSseq HDF5SummarizedExperiment"),
  make_option(c("-d","--in_dir"),type='character',default=".",
              help="Input file directory")
)
opt_parse=OptionParser(option_list=option_list)
opt = parse_args(opt_parse);
samples = unlist(strsplit(opt$input_samples,','));

#opt <- list()
#opt$output_bsseq_dir <- "testset_bs"
#opt$in_dir="~/work/projects/methylation-wdl/bismark_test"
#samples <- c("small_01", "small_02", "small_03")

message("Output sampleset directory: ", opt$output_bsseq_dir)
message("Input directory: ", opt$in_dir)
message("Samples: ", paste(samples, collapse=" "))

#########################
### Utility functions ###
#########################
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
				genome = "bisulfite genome of (.*) with",
				CHG_meth = "Total methylated C's in CHG context:\t(.*)",
				CHH_meth = "Total methylated C's in CHH context:\t(.*)",
				CHG_unmeth = "Total unmethylated C's in CHG context:\t(.*)",
				CHH_unmeth = "Total unmethylated C's in CHH context:\t(.*)"
				)

getPhenoData <- function(sample) {
  bismark_report <- file.path(opt$in_dir, paste0(sample, "_PE_report.txt"))
  report_lines <- readLines(bismark_report)
  pd <- data.frame(row.names=sample,
            t(as.data.frame(sapply(patterns, extract, report_lines), stringsAsFactors=FALSE)), stringsAsFactors = FALSE)
  # Convert columns with only numeric values from string to numeric
  suppressWarnings(convert_idx <- which(!is.na(sapply(colnames(pd), function(colnam) as.numeric(pd[,colnam])))))
  for (idx in convert_idx) pd[,idx] <- as.numeric(pd[,idx])
  return(pd)
}

getMethCov <- function(sample, gr) {
  tab <- read_tsv(file.path(opt$in_dir, paste0(sample, "_pe.bismark.cov.gz")),
                  col_types = "ciiiii",
                  col_names=c("chr", "pos", "pos2", "meth_percent", "m_count", "u_count"))
  tab_gr <- GRanges(tab$chr, IRanges(tab$pos, tab$pos))
  m <- rep(0, length(gr))
  cov <- rep(0, length(gr))
  ovl <- suppressWarnings(findOverlaps(tab_gr, gr))
  m[subjectHits(ovl)]   <- tab$m_count[queryHits(ovl)]
  cov[subjectHits(ovl)] <- tab$m_count[queryHits(ovl)] + tab$u_count[queryHits(ovl)]
  return(list(m=m, cov=as.integer(cov)))
}
#########################
#########################

# Set up genome-wide CpG GRanges
cpg <- DNAString("CG")
cpg_gr <- vmatchPattern(cpg, Hsapiens)
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

# Get phenodata
pd <- foreach(sample=samples, .combine="rbind") %do% getPhenoData(sample)

# Get methylation, coverage matrices
hdf5_m <- list()
hdf5_cov <- list()
#sample <- "small_02"
for(sample in samples) {
  message(sample)
  tmp <- getMethCov(sample, cpg_gr)
  m <- tmp$m
  cov <- tmp$cov
  hdf5_file <- paste0(sample, ".hdf5")
  if(file.exists(hdf5_file)) file.remove(hdf5_file)
  hdf5_m[[sample]] <- writeHDF5Array(matrix(m), name="m", file=hdf5_file)
  hdf5_cov[[sample]] <- writeHDF5Array(matrix(cov), name="cov", file=hdf5_file)
}

M <- do.call("cbind", hdf5_m)
Cov <- do.call("cbind", hdf5_cov)


# Save bsseq object to rds using HDF5Arrays for on-disk storage
bs <- BSseq(gr=cpg_gr, M=M, Cov=Cov, pData=pd, sampleNames=samples)
message("Saving HDF5SummarizedExperiment to ", opt$output_bsseq_dir)
#system.time(saveHDF5SummarizedExperiment(bs, dir=opt$output_bsseq_dir, chunk_dim=c(nrow(M), 1), replace = TRUE, verbose=TRUE))
system.time(saveHDF5SummarizedExperiment(bs, dir=opt$output_bsseq_dir, replace = TRUE, verbose=TRUE))
message("Done.")



