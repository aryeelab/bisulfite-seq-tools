library(optparse)
library(scmeth)


### UPDATE THIS TO CUSTOMIZE GENOME
library("BSgenome.Mmusculus.UCSC.mm10")
organism <- Mmusculus
genome <- "mm10"


option_list = list(
        make_option(c("-f","--bsFile"),action="store",type='character',default=NULL,
	                help="combined RDA file for entity set"),
	make_option(c("-o","--outDir"),action="store",type='character',default=NULL,
			help="output directory")
	)

opt_parse=OptionParser(option_list=option_list)
opt = parse_args(opt_parse)
bsfile <- opt$bsFile
outdir <- opt$outDir

load(bsfile)
report(combined_rda,outdir,organism,genome)