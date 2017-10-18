directory <- commandArgs(TRUE)[1]
bsgenome <- commandArgs(TRUE)[2]
library(dplyr)
library(scmeth)
library(bsgenome, character.only = TRUE)
bs <- SummarizedExperiment::loadHDF5SummarizedExperiment(directory)
report(bs, '.', get(bsgenome), unlist(lapply(strsplit(bsgenome, split="\\."), "[[",4)))
