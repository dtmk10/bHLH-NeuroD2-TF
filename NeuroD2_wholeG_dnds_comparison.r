library(dplyr)
library(ggplot2)

## 1) load and clean data
genome <- read.delim("~/documents/UPFClasses/PGB_Project/whole_genome_human_chimp.txt", sep="\t", header=TRUE)
colnames(genome) <- make.names(colnames(genome))

genome <- genome %>%
  mutate(omega = dN / dS) %>%
  filter(!is.na(omega) & omega < 10) %>%
  group_by(Ensembl.Gene.ID) %>%
  summarize(omega = mean(omega, na.rm=TRUE))

## at the end when we don't see NeuroD2 in the dataset, it's because it is filtered out here
## the reason it is filtered out is because it has a ds of 0
dnds <- dnds %>%
  mutate(omega = dN / dS) %>%
  filter(!is.na(omega) & omega < 10) 

### filtering here for unique values since the file we have has many duplicates. Next time download from ensemble
## with unique values only!
dnds_gene <- dnds %>%
  group_by(Ensembl.Gene.ID) %>%
  summarize(mean_dN = mean(dN, na.rm=TRUE),
            mean_dS = mean(dS, na.rm=TRUE),
            omega = mean(omega, na.rm=TRUE))


#### 2) filter data for bhlh factors only, and then combine into one dataset with filtered dndsvalue

bhlh_genes <- bhlh %>%
  filter(DBD=="bHLH") %>%
  select(Ensembl.ID)


#### build new dataset with unique values from dnds gene, and bhlh only genes.
#bhlh_dnds now holds the unique dnds values for the bhlh genes we want
bhlh_dnds <- dnds_gene %>%
  inner_join(bhlh_genes, by=c("Ensembl.Gene.ID"="Ensembl.ID"))

##select only neuroD2 for comparison
neurod2 <- bhlh_dnds %>%
  filter(Ensembl.Gene.ID == "ENSG00000171532")
neurod2 



## 3) Plot!
ggplot(bhlh_dnds, aes(x = omega)) +
  geom_histogram(bins = 40, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = neurod2$omega, color = "red", linewidth = 1.2) +
  annotate("text", x = neurod2$omega, y = Inf, label = "NeuroD2",
           vjust = 2, hjust = -0.1, color = "red") +
  labs(title = "dN/dS Ratios (ω) — Human vs Chimp bHLH Transcription Factors",
       x = "ω (dN/dS)", y = "Count") +
  theme_minimal(base_size = 12)




## diff type of plot
neuro <- as.numeric(neurod2$omega[1])
ggplot(bhlh_dnds, aes(x = omega)) +
  geom_histogram(
    aes(fill = after_stat(ifelse(xmin <= neuro & neuro < xmax, "NeuroD2 bin", "Other"))),
    bins = 40,
    color = "white",
    alpha = 0.8
  ) +
  geom_vline(xintercept = neuro, color = "red", linewidth = 1.2) +
  annotate("text", x = neuro, y = Inf, label = "NeuroD2", vjust = 1.5, color = "red") +
  scale_fill_manual(values = c("Other" = "steelblue", "NeuroD2 bin" = "red"), guide = "none") +
  labs(
    title = "dN/dS (ω) — bHLH TFs",
    x = "ω (dN/dS)", y = "Count"
  ) +
  theme_minimal(base_size = 12)

