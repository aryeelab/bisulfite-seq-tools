library(optparse)
library(scmeth)


option_list = list(
        make_option(c("-f","--bsFile"),action="store",type='character',default=NULL,
	                help="combined RDA file for entity set"),
	make_option(c("-o","--outDir"),action="store",type='character',default=NULL,
			help="output directory"),
	#make_option(c('-r','--organism'),action="store",type='character',default=NULL,
	#		help='name of organism'),
	#make_option(c('-g','--genome'),action="store",type='character',default=NULL,
	#		help='reference genome')
	)

opt_parse=OptionParser(option_list=option_list)
opt = parse_args(opt_parse)
bsfile <- opt$bsFile
outdir <- opt$outDir
#organism <- opt$organism
#genome <- opt$genome


load(bsfile)


### UPDATE THIS TO CUSTOMIZE GENOME
library("BSgenome.Mmusculus.UCSC.mm10")
organism <- Mmusculus
genome <- "mm10"


report(combined_rda,outdir,organism,genome)