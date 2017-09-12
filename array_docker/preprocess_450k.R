
library(optparse)
option_list = list(
        make_option(c("-d","--data_dir"),action="store",type='character',default=NULL,
	                help="directory where data files are"),
	make_option(c("-s","--sample_sheet"),action="store",type='character',default=NULL,
			help="CSV describing sample names, attributes and files"),
	make_option(c('-q','--qc_intensity_cutoff'),action="store",type='character',default=NULL,
			help='QC intensity cutoff'),
        make_option(c('-b','--beta_file'),action="store",type='character',default=NULL,
	                help='output file with beta values'),
	make_option(c('-c','--qc_file'),action="store",type='character',default=NULL,
			help='output file for QC data'),
        make_option(c('-p','--preprocessed_genomic_ratio_set_rds'),action="store",type='character',default=NULL,
	                help='preprocessed genomic ratio set RDS')
	)
opt_parse=OptionParser(option_list=option_list)
opt = parse_args(opt_parse)

library(minfi)
# INPUTS
cat('started processing inputs...')
idat_base_dir <- opt$data_dir  #<- "testdata/450k"
sample_sheet  <- opt$sample_sheet #<- "testdata/450k/samples.csv"
qc_intensity_cutoff <- as.numeric(opt$qc_intensity_cutoff)  #<- 11
cat('finished processing inputs\n')

# OUTPUT
cat('started processing outputs...')
qc_file <- opt$qc_file #<- "qc.tsv"
betas_file <- opt$beta_file #<- "beta.tsv"
preprocesed_genomic_ratio_set_rds <- opt$preprocessed_genomic_ratio_rds #<- "gmset.rds"
cat('finished processing outputs\n')

cat('reading targets\n')
targets <- read.csv(sample_sheet, stringsAsFactors = FALSE)
targets$Basename <- file.path(idat_base_dir, targets$Basename)

# This reads the files specified in the Basename column of targets
rgset <- read.metharray.exp(targets=targets)
sampleNames(rgset) <- targets[,1]

# QC
cat('starting QC...')
mset <- preprocessRaw(rgset) 
qc <- cbind(pData(rgset), getQC(mset))
bad_idx <- qc$mMed < qc_intensity_cutoff
badsamples <- rownames(qc)[bad_idx] 
qc$qc_pass <- ifelse(sampleNames(mset) %in% badsamples, "FAIL", "PASS")
plotQC(qc, badSampleCutoff = qc_intensity_cutoff)
cat('Done with QC\n')

# Preprocess and get betas
grset_funnorm <- preprocessFunnorm(rgset[,!bad_idx], nPCs=2)  
b <- getBeta(grset_funnorm)

# Save output files
write.table(qc, file=qc_file, sep="\t", quote=FALSE, row.names=FALSE)
write.table(b, file=betas_file, sep="\t", quote=FALSE, col.names=NA)
saveRDS(grset_funnorm, file=preprocesed_genomic_ratio_set_rds)  

