
# This script combines multiple bisulfite sequencing samples processed by the Bismark aligner
# into a single BSseq object saved as an HDF5SummarizedExperiment.
#

required_libraries <- c("readr","stringr","HDF5Array","bsseq","foreach","optparse","Biostrings","GenomeInfoDb")
sapply(required_libraries, library, character.only=TRUE)

options(DelayedArray.block.size=4500e6)



option_list = list(
   make_option(c("-i","--input_pe_report_files"),type='character',default=NULL,
               help="Output from bismark_methylation_extractor. Use comma-separated file names "),
   make_option(c("-j","--input_covgz_files"),type='character',default=NULL,
               help="Output from bismark_methylation_extractor. Use comma-separated sample names "),
   make_option(c("-k","--input_mbias_files"),type='character',default=NULL,
               help="Output from bismark_methylation_extractor. Use comma-separated file names "),
   make_option(c("-o","--output_bsseq_dir"),type='character',default=NULL,
               help="Directory in which to save the BSseq HDF5SummarizedExperiment"),
   make_option(c("-m","--mbias_dir"),type='character',default=NULL,
               help="Directory in which to save the mbias files"),
   make_option(c("-d","--in_dir"),type='character',default=".",
               help="Input file directory"),
   make_option(c("-g","--bsgenome"),type='character',default=".",
               help="Genome package")
 )
opt_parse=OptionParser(option_list=option_list)
opt = parse_args(opt_parse);
library(opt$bsgenome, character.only=TRUE)
pe_report_files = scan(opt$input_pe_report_files, what="character")
covgz_files = scan(opt$input_covgz_files,what="character")
mbias_files = scan(opt$input_mbias_files,what="character")



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

getPhenoData <- function(pe_report_file) {
  bismark_report <- file.path(pe_report_file)
  samplename <- basename(bismark_report)
  samplename <- gsub(samplename, pattern="_PE_report.txt", replacement="")
  report_lines <- readLines(bismark_report)
  pd <- data.frame(row.names=samplename,
                   t(as.data.frame(sapply(patterns, extract, report_lines), stringsAsFactors=FALSE)), stringsAsFactors = FALSE)
  # Convert columns with only numeric values from string to numeric
  suppressWarnings(convert_idx <- which(!is.na(sapply(colnames(pd), function(colnam) as.numeric(pd[,colnam])))))
  for (idx in convert_idx) pd[,idx] <- as.numeric(pd[,idx])
  return(pd)
}

getMethCov <- function(covgz_file, gr) {
  tab <- read_tsv(covgz_file,
                  col_types = "ciidii",
                  col_names=c("chr", "pos", "pos2", "meth_percent", "m_count", "u_count"))
  tab_gr <- GRanges(tab$chr, IRanges(tab$pos, tab$pos))
  m <- rep(0, length(gr))
  cov <- rep(0, length(gr))
  ovl <- suppressWarnings(findOverlaps(tab_gr, gr))
  m[subjectHits(ovl)]   <- tab$m_count[queryHits(ovl)]
  cov[subjectHits(ovl)] <- tab$m_count[queryHits(ovl)] + tab$u_count[queryHits(ovl)]
  return(list(m=m, cov=as.integer(cov)))
}



transferMbias <- function(mbias_file) {
  fileName<-basename(mbias_file)
  newPath<-paste0(opt$mbias_dir,'/',fileName)
  file.copy(from=mbias_file,to=newPath)
}


#########################
#########################

cpg_gr <- DNAString("CG")
cpg_gr <- vmatchPattern(cpg_gr, get(opt$bsgenome))
cpg_gr <- keepStandardChromosomes(cpg_gr, pruning.mode="coarse")


# Set up genome-wide CpG GRanges
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
cat(pe_report_files,sep="\n")
pd <- foreach(pe_report_file=pe_report_files, .combine="rbind") %do% getPhenoData(pe_report_file)
rownames(pd) <- sub("_report.txt", "", rownames(pd))

# Store the mbias file in a given directory
cat(mbias_files, sep='\n')
foreach(mbias_file=mbias_files) %do% transferMbias(mbias_file)


# Get methylation, coverage matrices
hdf5_m <- list()
hdf5_cov <- list()
for(covgz_file in covgz_files) {
  tmp <- getMethCov(covgz_file, cpg_gr)
  m <- tmp$m
  cov <- tmp$cov
  samplename <- basename(covgz_file)
  samplename <- gsub(samplename, pattern="_PE_report.txt", replacement="")
  hdf5_file <- paste0(samplename, ".hdf5")
  if(file.exists(hdf5_file)) file.remove(hdf5_file)
  hdf5_m[[samplename]] <- writeHDF5Array(matrix(m), name="m", file=hdf5_file)
  hdf5_cov[[samplename]] <- writeHDF5Array(matrix(cov), name="cov", file=hdf5_file)
}
message("Generated coverage and methylation info")

M <- do.call("cbind", hdf5_m)
Cov <- do.call("cbind", hdf5_cov)


# Save bsseq object to rds using HDF5Arrays for on-disk storage
bs <- BSseq(gr=cpg_gr, M=M, Cov=Cov, pData=pd, sampleNames=rownames(pd))
message("Saving HDF5SummarizedExperiment to ", opt$output_bsseq_dir)
#system.time(saveHDF5SummarizedExperiment(bs, dir=opt$output_bsseq_dir, chunk_dim=c(nrow(M), 1), replace = TRUE, verbose=TRUE))
system.time(saveHDF5SummarizedExperiment(bs, dir=opt$output_bsseq_dir, replace = TRUE, verbose=TRUE))
message("Done.")
