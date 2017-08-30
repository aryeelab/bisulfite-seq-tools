library(minfi)

# INPUTS
idat_base_dir <- "testdata/450k"
sample_sheet <- "testdata/450k/samples.csv"
qc_intensity_cutoff <- 11

# OUTPUT
qc_file <- "qc.tsv"
betas_file <- "beta.tsv"
preprocesed_genomic_ratio_set_rds <- "gmset.rds"

targets <- read.csv(sample_sheet, stringsAsFactors = FALSE)
targets$Basename <- file.path(idat_base_dir, targets$Basename)

# This reads the files specified in the Basename column of targets
rgset <- read.metharray.exp(targets=targets)
sampleNames(rgset) <- targets[,1]

# QC
mset <- preprocessRaw(rgset) 
qc <- cbind(pData(rgset), getQC(mset))
bad_idx <- qc$mMed < qc_intensity_cutoff
badsamples <- rownames(qc)[bad_idx] 
qc$qc_pass <- ifelse(sampleNames(mset) %in% badsamples, "FAIL", "PASS")
plotQC(qc, badSampleCutoff = qc_intensity_cutoff)

# Preprocess and get betas
grset_funnorm <- preprocessFunnorm(rgset[,!bad_idx], nPCs=2)  
b <- getBeta(grset_funnorm)

# Save output files
write.table(qc, file=qc_file, sep="\t", quote=FALSE, row.names=FALSE)
write.table(b, file=betas_file, sep="\t", quote=FALSE, col.names=NA)
saveRDS(grset_funnorm, preprocesed_genomic_ratio_set_rds)  

