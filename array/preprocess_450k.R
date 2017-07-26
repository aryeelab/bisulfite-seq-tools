library(minfi)

# INPUTS
idat_base_dir <- "testdata/450k"
sample_sheet <- "testdata/450k/samples.csv"
qc_intensity_cutoff <- 11

# OUTPUT
qc_file <- "qc.csv"
preprocesed_genomic_ratio_set <- "gmset.rds"

targets <- read.csv(sample_sheet, stringsAsFactors = FALSE)
targets$Basename <- file.path(idat_base_dir, targets$Basename)

# This reads the files specified in the Basename column of targets
rgset <- read.metharray.exp(targets=targets)

# QC
mset <- preprocessRaw(rgset) 
qc <- cbind(pData(rgset), getQC(mset))
bad_idx <- qc$mMed < qc_intensity_cutoff
badsamples <- rownames(qc)[bad_idx] 
qc$qc_pass <- ifelse(sampleNames(mset) %in% badsamples, "FAIL", "PASS")
plotQC(qc, badSampleCutoff = cutoff)
write.csv(qc, file=qc_file)

# Preprocess
grset_funnorm <- preprocessFunnorm(rgset[,!bad_idx], nPCs=2)  
saveRDS(grset_funnorm, preprocesed_genomic_ratio_set)  

