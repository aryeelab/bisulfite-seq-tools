library("optparse")
library("scmeth")
library("bsseq")
library(rmarkdown)

### UPDATE THIS TO CUSTOMIZE GENOME
library("BSgenome.Mmusculus.UCSC.mm10")
organism <- Mmusculus
genome <- "mm10"


option_list = list(
        make_option(c("-i","--inputBsseqHDF5"),action="store",type='character',default=NULL,
	                help="Directory with bsseq HDF5SummarizedExperiment"),
	make_option(c("-o","--outDir"),action="store",type='character',default=NULL,
			help="output directory")
	)

opt_parse=OptionParser(option_list=option_list)
opt = parse_args(opt_parse)

bs <- loadHDF5SummarizedExperiment(opt$inputBsseqHDF5)
report(bs, opt$outDir, organism, genome)
