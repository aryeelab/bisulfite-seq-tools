library(ggplot2)
df <- read.table("firecloud_runtimes.txt", header=TRUE)
colnames(df)[1] <- "VM"
head(df)
p <- ggplot(df, aes(hours_per_1B_reads_30X, cost_per_1B_reads_30X)) + geom_point(aes(color=VM)) + xlim(0, NA) + ylim(0, NA) + theme_bw() + xlab("Processing time (hours)") + ylab("Cost ($)")
p
ggsave("firecloud_runtimes.pdf", p, width = 6, height=3)
