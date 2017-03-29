library(bsseq)
library(Biobase)
cpgtab <- Sys.getenv("CPGTAB")
rda <- Sys.getenv("RDA")

message("Reading methylation evidence from ", cpgtab)
pd <- data.frame(row.names=sub(".bedGraph", "", basename(cpgtab)), type=sub("-.*.bedGraph", "", basename(cpgtab)), stringsAsFactors=F)
print(head(pd))


#tab <- read.delim(cpgtab, nrows=100000)
tab <- read.delim(cpgtab,skip=1)
print(head(tab))
m <- tab[, 5] # methylated count
um <- tab[, 6] # unmethylated count
cov <- m+um # total coverage
idx <- m>cov
if (sum(idx)>0) {
    message("Warning: ", sum(idx), " C count(s) are greater than the CT count")
    m[idx] <- cov[idx]
}
bs <- BSseq(chr=tab[,1], pos=tab[,3], M=matrix(m), Cov=matrix(cov), pData=pd, sampleNames=pd$sample_id)

message("Saving bs to ", rda)
save(bs, file=rda)
message("Done.")

