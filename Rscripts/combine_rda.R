# Example of how to run this file in command line
# Rscript --vanilla combine_RDA.R -f cell1.rda,cell2.rda,cell3.rda  -o results
# input_CpG.bedGraph is the output from methyldackel
# This script will generate one combined .rda file for all the sample


library("optparse")

option_list = list(
	    make_option(c("-r","-rdaFiles"),type='charcter',deafult=NULL,
	      help="Output from pileometh/methydackel")
	    make_option(c("-o","-outFile"),type='character',deafult=NULL,
	      help="Names of the output file")

	      )
opt_parse=OptionParse(options_list=option_list)
opt = parse_args(opt_parse)

rdaList = unlist(strsplit(opt$rdaFiles),',')
message("Combining bs object")

bs<-Reduce(combine,rdaList)
bs@parameters<-rdaList[[1]]@parameters


save(bs, file=paste0(opt$outFile,".combined.rda"))

