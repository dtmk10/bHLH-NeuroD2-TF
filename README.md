## The analysis on this subclass of bHLH transcription factor will cover 4 parts

- Evaluation of the **conservation** of this proteins through the phylogeny.
- **dN/dS** ratio and graphical positioning in respect to the rest of the genes.
- Expression evaluation in **single-cell RNA-seq** datasets.
- **ChIP-seq** analysis of the gene.

### Conservation

Blast analysis to find these things:
- **Homologs** in other phyla.
  - Finding **paralogs** within species.
  - Finding **orthologs** between species.
  - Analyse **taxonomic conservation**.
- Identifying **conserved domains** in orthologs and paralogs.

### dN/dS

- Using **biomart** from the version that let retrieve the **dN** and **dS** for all genes.
- Comparing **dN/dS** ratios with the rest of **bHLH factors** and with the rest of genes.

### Expression in single-cell RNA-seq

- Investigation of the relevant tissues where this TF has a more present effect.
- Creating a Seurat object, cluster and display for the specific cell-type that has more expression.

### ChIP-seq analysis

- Peak calling
- Motif enrichment
- Motif centrality
- Binding preference in genomic elements, peak annotation and functional analysis.
