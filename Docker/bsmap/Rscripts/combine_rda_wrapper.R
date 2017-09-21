# Example of how to run this file in command line
# Rscript --vanilla combine_rda_wrapper.R -f cell1.rda,cell2.rda,cell3.rda  -o results
# input_CpG.bedGraph is the output from methyldackel
# This script will generate one combined .rda file for all the sample


library("optparse")
library(BiocGenerics)
library(scmeth)
option_list = list(
            make_option(c("-r","--rdaFiles"),type='character',default=NULL,
              help="Output from pileometh/methydackel"),
            make_option(c("-o","--outFile"),type='character',default=NULL,
              help="Names of the output file")
              )
opt_parse=OptionParser(option_list=option_list)
opt = parse_args(opt_parse);
rdaList = unlist(strsplit(opt$rdaFiles,','));
message("Combining bs object");
#bslist <- vector(mode="character", length=length(rdaList))
for(i in 1:length(rdaList)){
      load(rdaList[i])
      if(!exists('bslist')){
            bslist <- bs
      }else{
        bslist = c(bslist,bs);
      }
}
combined_rda <- combineRDA(bslist)
save(combined_rda, file=paste0(opt$outFile,".combined.rda"))




