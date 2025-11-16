# The analysis on this subclass of bHLH transcription factor will cover 4 parts

- **BLAST** analysis for conservation across species.
- **dN/dS** ratio and graphical positioning in respect to the rest of the genes.
- **Single-cell RNA-seq** expression analysis of the gene.
- **ChIP-seq** analysis of the gene.
# Folder structure
- `fasta_files/`: contains the fasta files with the sequences used in the different analysis.
- `Blast_Conservation/`: contains the scripts and results for the BLAST conservation analysis.
- `dN_dS/`: contains the scripts and results for the dN/dS analysis.
- `scRNA_seq/`: contains the scripts and results for the single-cell RNA-seq
- `ChIP_seq/`: contains the scripts and results for the ChIP-seq analysis.
- `figures/`: contains all the figures generated during the analysis.
# Analysis steps
### BLAST conservation
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



# Good Luck 
