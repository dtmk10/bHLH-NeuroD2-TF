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


gene <- seq2gene(peaksResized, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
pathway2 <- enrichPathway(gene, organism = "mouse")

length(pathway2@gene)
sum(pathway2@result$p.adjust < 0.05)
gene_list_adjusted <<- "x"
for(i in 1:length(pathway2@gene)){
  if(pathway2@result$p.adjust[i] < 0.05){
    gene_list_adjusted <<- paste(gene_list_adjusted, pathway2@gene[i], sep = "\n")
  }
}
write(gene_list_adjusted, file = "geneTable.txt", append = FALSE)

dotplot(pathway2)
barplot(pathway2)

