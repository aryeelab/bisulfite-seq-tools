# Example of how to run this file in command line
# Rscript --vanilla create_RDA.R -f input_CpG.bedGraph -o results
# input_CpG.bedGraph is the output from methyldackel
# This script will generate the .rda for that file

library(bsseq)
library(Biobase)
library("optparse")

option_list = list(
	    make_option(c("-f","--cpgFile"),type='character',default=NULL,
	      help="Output from pileometh/methydackel"),
	    make_option(c("-o","--outFile"),type='character',default=NULL,
	      help="Names of the output file")

	      )
opt_parse=OptionParser(option_list=option_list)
opt = parse_args(opt_parse)

tab<-read.delim(opt$cpgFile,sep='\t',header=F);
tab<-tab[!is.na(tab$V2),];
message('Creating BSseq object with ', nrow(tab),' loci.')
sample<-sub("_CpG.bedGraph","",basename(opt$cpgFile))
m<-tab[,5] # methylated count
um<-tab[,6] # unmethylated count
cov<-m+um

message("generating bs object")
bs<-BSseq(chr=tab[,1], pos=tab[,3], M= matrix(m), Cov=matrix(cov),sampleNames=sample)
save(bs, file=paste0(opt$outFile,".rda"))
message("Done.")


