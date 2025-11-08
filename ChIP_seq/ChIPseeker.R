library("ChIPseeker")
library("ReactomePA")
# annotation of Mus musculus
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
# annotation of Mus musculus in ENTREZ format
# Used in the functional analysis
library(org.Mm.eg.db)
mmdb <- org.Mm.eg.db

browseVignettes("ChIPseeker")

setwd("~/Documents/PGB-Workspace/bHLH-project/ChIP_seq/PeakCalling")

# The input is the peaks.narrowPeaks from macs2
narrowPeaks <- readPeakFile("IP_peaks.narrowPeak")
# The resized peaks with the middle point adjusted
peaksResized <- readPeakFile("IP_peaks.resized.bed")


# Peak spread over the genome
covplot(narrowPeaks, weightCol = "V5")

# Specifying chromosome and coords to visualize
covplot(narrowPeaks, weightCol="V5", chrs=c("chr1", "chr18"), xlim=c(4.5e7, 5e7))


## ==========================
# ChIP-PROFILING
## ==========================

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

tagMatrix <- getTagMatrix(peaksResized, windows=promoter)

peak_Profile_Heatmap(peak = peaksResized,
                     upstream = 1000,
                     downstream = 1000,
                     by = "gene",
                     type = "start_site",
                     TxDb = txdb,
                     nbin = 800)



## ===================
# PEAK ANNOTATION
## ===================

peakAnno <- annotatePeak(peaksResized, tssRegion=c(-3000, 3000),
                         TxDb=txdb)

plotAnnoPie(peakAnno)
vennpie(peakAnno)
plotDistToTSS(peakAnno, title = "Dist to TSS")




## ==================
# FUNCTIONAL ENRICHMENT ANALYSIS using reactomepa
## ==================

library(org.Mm.eg.db)
mmdb <- org.Mm.eg.db


pathway1 <- enrichPathway(as.character(as.data.frame(peakAnno)$geneId), organism = "mouse")

sum(pathway1@result$p.adjust < 0.05)

dotplot(pathway1)


