

setwd("~/Documents/PGB-Workspace/bHLH-project/ChIP_seq/IP_meme/")

memeSummary <- read.table(file = "summary.tsv", header = TRUE)
centrimoTsv <- read.table(file = "centrimo_out/centrimo_out_fixed.tsv", sep = "\t", header = TRUE)

