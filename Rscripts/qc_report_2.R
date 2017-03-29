



knitc = function(input) {
  library(knitr)
  opts_chunk$set(cache.path = paste('cache/QC', input, sep = ''))
  knit(input)
}
library(ggplot2)
library(reshape2)
library(bsseq)
library(matrixStats)
library(AnnotationHub)
library(BiocGenerics,lib.loc='/apps/lab/aryee/R/R-3.2.3/lib64/R/site-library')
library(GenomicRanges,lib.loc='/apps/lab/aryee/R/R-3.2.3/lib64/R/site-library')
library(foreach)
library(scales)
library(Repitools)
library(BSgenome.Hsapiens.UCSC.hg19,lib.loc='/apps/lab/aryee/R/R-3.2.3/lib64/R/site-library')
library(data.table)
#library(Matrix)
if (basename(getwd())!="code") setwd("code")

metrics<-Sys.getenv("MET")
bsseq<-Sys.getenv("BS")








#' ## Read and mapping metrics

#+ echo=FALSE, fig.width=8, fig.height=7,cache.lazy=FALSE 
load(paste0(Dir,bsseq))




#+ echo=FALSE
#dat <- read.delim(metrics, header = TRUE,sep=" ", row.names=1)
# Plot mapped and unmapped reads
dat <- read.delim(paste0(Dir,metrics), header = TRUE,sep=" ", row.names=1)
dat$sample <- basename(rownames(dat))
dat$Unmapped<-dat$total - dat$mapped
CellDF<-DF<-data.frame(matrix(unlist(strsplit(dat$sample,"[-]")),nrow = length(strsplit(dat$sample,"[-]")),byrow=T))
dat$class<-CellDF[,1]
m <- melt(dat[,c("sample", "mapped", "Unmapped")], id.vars = "sample", variable.name = "Mapping_status")
o <- order(dat$class, dat$total)
sampleOrder <- dat$sample[o]
m$sample <- factor(m$sample, levels=sampleOrder)
Samples<-grep('bulk',m$sample,invert=TRUE)
ggplot(m[Samples,], aes(sample, value, fill=Mapping_status)) + geom_bar(stat="identity") + coord_flip() + scale_y_continuous(name="Number of reads", labels = comma) + xlab("Sample") + ggtitle("Number of reads in Samples") + theme_bw()

#' ## CpG Coverage
#+ echo=FALSE,cache=TRUE, cache.lazy=FALSE
#load(bsseq)
seqlevels(bs)<-paste0("chr",seqlevels(bs))
keep <- seqlevels(bs) %in% names(Hsapiens)
#seqlevels(bs) <- seqlevels(bs)[keep]
idx <- seqnames(bs) %in% seqlevels(bs)[keep]
bs <- bs[idx,]
seqlevels(bs)<-as.character(unique(seqnames(bs)))
cov <- getCoverage(bs)
Samples<-grep("bulk",colnames(cov),invert=TRUE)
cov<-cov[,Samples]
CellDF<-DF<-data.frame(matrix(unlist(strsplit(colnames(cov),"[-]")),nrow = length(strsplit(colnames(cov),"[-]")),byrow=T))
DFSamples<-grep("bulk",CellDF[,2],invert=TRUE)
CellDF<-CellDF[DFSamples,]
sampleName<-colnames(cov)
covDf <- data.frame(sample=sampleName, coveredCpgs=colSums(cov>=1))
covDf$grp <- substr(covDf$sample, 1, 1)
AllCpG_Assigned<- melt(covDf, id.vars = c("grp", "sample"))
m <- melt(covDf, id.vars = c("grp", "sample"))
o <- order(covDf$grp, covDf$coveredCpgs)
sampleOrder <- covDf$sample[o]
m$sample <- factor(m$sample, levels=sampleOrder)
ggplot(m[Samples,], aes(sample, value)) + geom_hline(data = m, aes(yintercept = as.numeric(100000)), linetype=2, color="blue")+ geom_bar(stat="identity") + coord_flip() + scale_y_continuous(name="Number of CpGs", labels = comma) + xlab("Sample") + ggtitle("Number of Covered CpGs")  + theme_bw()



