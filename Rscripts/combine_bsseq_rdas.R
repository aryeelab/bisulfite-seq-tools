library(bsseq)
library(Biobase)

infiles <- Sys.getenv("INFILES")
infiles <- unlist(strsplit(infiles, " "))
outfile <- Sys.getenv("OUTFILE")

message("Reading bsseq objects from:\n", paste(infiles, collapse="\n"))
bsList <- lapply(infiles, function(file) {
    load(file)
    bs
})

bs <- Reduce(combine, bsList)

# Record smoothing parameters (if any)
bs@parameters <- bsList[[1]]@parameters

# Add colors for plotting
pData <- pData(bs)
pData$col <- as.numeric(factor(pData$type))
pData(bs) <- pData

message("Saving bs to ", outfile)
save(bs, file=outfile)
message("Done.")

