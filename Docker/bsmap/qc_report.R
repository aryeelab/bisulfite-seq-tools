#+ echo=FALSE, message=FALSE
sampleset <- "2017-04-21-SN0118832"
bsFile <- "/data/aryee/scrrbs/flowcells/2017-04-21_SN0118832/2017-04-21-SN0118832.bsseq.rda"
metrics_file <- "/data/aryee/scrrbs/flowcells/2017-04-21_SN0118832/2017-04-21-SN0118832.read_metrics.txt"

#+ echo=FALSE, message=FALSE
pipeline_ver <- Sys.getenv("PIPELINE_VER", "v0")
genome<-tolower(Sys.getenv("GENOME"))
genome <- "grch38"
valid_barcodes_file <- "/data/aryee/ma695/work/projects/scdnam/valid_barcodes_scrrbs.txt"

#+ echo=FALSE, message=FALSE
library(yaml)
library(ggplot2)
library(reshape2)
library(bsseq)
library(matrixStats)
library(foreach)
library(scales)
library(Repitools)
library(data.table)
library(stringr)
library(knitr)
#opts_chunk$set(fig.path=file.path(libdir, "reports/qc_report_figures/"))
valid_barcodes <- read.table(valid_barcodes_file, stringsAsFactors = FALSE)[,1]

#+ load_genome_name_from_config, echo=FALSE, message=FALSE
config <- yaml.load_file("../config.yaml")
bsPkgName <- config$genomes[[genome]]$bsgenome
library(bsPkgName, character.only = TRUE)
bsgenome <- eval(parse(text = strsplit(bsPkgName, split="\\.")[[1]][2]))

#' ## Number of reads
#+ number_of_reads, echo=FALSE, fig.width=9, fig.height=5, dpi=300, fig.retina=2
dat <- read.delim(metrics_file, header = TRUE,sep=" ", row.names=1)
#dat$sample <- sub(".*_(.*)", "\\1", dat$sample)
dat$sample <- rownames(dat)
dat$sample_short <- str_split_fixed(dat$sample, "_", 2)[,2]
dat$pool_id <- str_split_fixed(dat$sample, "_", n=2)[,1]
#x<-sapply(valid_barcodes, grepl, dat$sample)
#rowSums(sapply(valid_barcodes, grepl, dat$sample))
dat$Unmapped<-dat$total - dat$mapped
m <- melt(dat[,c("sample", "sample_short", "pool_id", "mapped", "Unmapped")], id.vars = c("sample", "sample_short", "pool_id"), variable.name = "Mapping_status")
o <- order(dat$total)
sampleOrder <- dat$sample[o]
m$sample <- factor(m$sample, levels=sampleOrder)
ggplot(m, aes(sample, value/1e6, fill=Mapping_status)) + geom_bar(stat="identity") + coord_flip() + scale_y_continuous(name="Number of reads (millions)", labels = comma) + xlab("Sample") + ggtitle(paste0(sampleset, "\nNumber of reads")) + theme_bw(base_size = 12)

#+ number_of_reads_by_pool, echo=FALSE, fig.width=9, fig.height=5, dpi=300, fig.retina=2
m$sample_short <- factor(m$sample_short, levels=names(sort(tapply(dat$total, dat$sample_short, sum))))
ggplot(m, aes(sample_short, value/1e6, fill=Mapping_status)) + geom_bar(stat="identity") + coord_flip() + facet_wrap(~pool_id, scales="free_y") + scale_y_continuous(name="Number of reads (millions)", labels = comma) + xlab("Sample") + ggtitle(paste0(sampleset, "\nNumber of reads")) + theme_bw(base_size = 12)

#' ## CpG Coverage
#+ covered_cpgs, echo=FALSE, fig.width=9, fig.height=5, dpi=300, fig.retina=2
#bsFile <- file.path(dir, paste0(sampleset, ".bsseq.rda"))
#bsFile <- "/PHShome/ma695/work/projects/sc_dnmt_ko/output/2017-04-21_gp_rrbs.rda"
load(bsFile)
# Temporarily removing this line until we fix the genome slot of the gr in the bs object (in create_bsseq_rda.R)
#bs <- keepStandardChromosomes(bs)
cov <- getCoverage(bs)
#sampleName <- sub(".*_(.*)", "\\1", colnames(cov))
covDf <- data.frame(sample=colnames(cov), 
                    pool_id <- str_split_fixed(colnames(cov), "_", n=2)[,1],
                    sample_short=str_split_fixed(colnames(cov), "_", n=2)[,2],
                    coveredCpgs=colSums(cov>=1))
o <- order(covDf$coveredCpgs)
sampleOrder <- covDf$sample[o]
covDf$sample <- factor(covDf$sample, levels=sampleOrder)
ggplot(covDf, aes(sample, coveredCpgs)) + geom_hline(data = m, aes(yintercept = as.numeric(1000000)), linetype=2, color="blue")+ geom_bar(stat="identity") + coord_flip() + scale_y_continuous(name="Number of CpGs", labels = comma) + xlab("Sample") + ggtitle(paste0(sampleset, "\nNumber of Covered CpGs"))  + theme_bw(base_size = 12)

#+ covered_cpgs_by_pool, echo=FALSE, fig.width=9, fig.height=5, dpi=300, fig.retina=2
ggplot(covDf, aes(pool_id, coveredCpgs)) + geom_hline(data = m, aes(yintercept = as.numeric(1000000)), linetype=2, color="blue")+ geom_boxplot() + coord_flip() + scale_y_continuous(name="Number of CpGs", labels = comma) + xlab("Sample") + ggtitle(paste0(sampleset, "\nNumber of Covered CpGs (by pool)"))  + theme_bw(base_size = 12)
#table(covDf$coveredCpgs>=100000)
#x <- subset(covDf, coveredCpgs>=100000)
#median(x$coveredCpgs)

#' ### Median number of covered CpGs
median(covDf$coveredCpgs)
tapply(covDf$coveredCpgs, covDf$sample_short, median)

#' ## Sequencing saturation analysis
#+ saturation_analysis, eval=TRUE, echo=FALSE, fig.width=7, fig.height=5, dpi=300, fig.retina=2
dsRates <- seq(0.01, 1, 0.1)
dsCov <- apply(cov, 2, function(x) {
  t <- table(x)
  n <- as.numeric(names(t))
  foreach(dsRate=dsRates, .combine="c") %do% {
    obsRate <- 1-pbinom(0, n, dsRate)
    head(t)
    sum(t[n>0])
    dsT <- t*obsRate
    head(dsT)
    dsCoveredCpgs <- sum(dsT[n>0])
    dsCoveredCpgs
  }
})
rownames(dsCov) <- dsRates
#colnames(dsCov) <- sub(".*_(.*)", "\\1", colnames(dsCov))
m <- melt(dsCov, varnames=c("downsampleRate", "cell"), value.name="coverage")
m$reads <- m$downsampleRate * sum(dat$total)
ggplot(m, aes(reads/1e6, coverage, group=cell)) + geom_line() + theme_bw() + xlab("Total library pool reads (millions)") + ylab("CpG coverage") + ggtitle(paste0(sampleset, "\nSequencing saturation analysis"))



#+ echo=FALSE
if (FALSE) {
  ###################################################################################
  ## Repeat-masked coverage dropped until we implement configurable genome support ##
  ###################################################################################
  #' CpGs in non-repeat masker regions
  #+ covered_nonrepeatmasked_cpgs, echo=FALSE, fig.width=9, fig.height=1+sum(covLogicVector)/5, fig.retina=2
  #seqlevels(repeatGr) <- sub("chr", "", seqlevels(repeatGr)) # need to remove chr in order to match
  rep <- countOverlaps(bs, repeatGr)>0
  cov_Nrep <- cov[!rep,covLogicVector]
  sampleName<-colnames(cov_Nrep)
  covDf <- data.frame(sample=sampleName, coveredCpgs=colSums(cov_Nrep>=1))
  #AllCpG <- melt(covDf, id.vars = "sample")
  m <- melt(covDf, id.vars = "sample")
  #Samples<-grep('bulk',m$sample,invert=TRUE)
  o <- order(covDf$coveredCpgs)
  sampleOrder <- covDf$sample[o]
  m$sample <- factor(m$sample, levels=sampleOrder)
  ggplot(m, aes(sample, value)) +  geom_bar(stat="identity") + coord_flip() + scale_y_continuous(name="Number of CpGs", labels = comma) + xlab("Sample") + ggtitle("Number of Covered Non-repeatmasked CpGs")  + theme_bw()
}