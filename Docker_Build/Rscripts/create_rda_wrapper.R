library(bsseq)
library(Biobase)
library(optparse)
library(scmeth)
option_list = list(
        make_option(c("-f","--cpgFile"),action="store",type='character',default=NULL,
	        help="Output from pileometh/methydackel"),
	make_option(c("-o","--outFile"),action="store",type='character',default=NULL,
		help="Names of the output file"),
	make_option(c('-r','--read.metrics.file'),action="store",type='character',default=NULL,
		help='read metrics file'),
	make_option(c('-b','--bsconversion'),action="store",type='character',default=NULL,
		help='bisulphite conversion rate')
	)

opt_parse=OptionParser(option_list=option_list)
opt = parse_args(opt_parse)

cpgfile <- opt$cpgFile
readmetrics <- opt$read.metrics.file
outfile <- opt$outFile
bsconversionfile <- opt$bsconversion

bs <- createRDA(cpgfile,readmetrics,bsconversionfile);

save(bs,file=outfile)