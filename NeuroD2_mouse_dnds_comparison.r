library(dplyr)
library(ggplot2)

## 1) Load and clean data
mouse_dnds <- read.delim("/Users/austingilbride/PGB/bHLH-NeuroD2-TF/NeuroD2_bhlh_human_mouse.txt")
bhlh <- read.delim("~/documents/UPFClasses/PGB_Project/HumanTFs_DBD.txt", sep="\t", header=TRUE)
genome <- read.delim("/Users/austingilbride/PGB/bHLH-NeuroD2-TF/human_mouse_whole_genome.txt", sep="\t", header=TRUE)
colnames(genome) <- make.names(colnames(genome))


genome <- genome %>%
  mutate(omega = dN / dS) %>%
  filter(!is.na(omega) & omega < 10) %>%
  group_by(Ensembl.Gene.ID) %>%
  summarize(omega = mean(omega, na.rm=TRUE))


mouse_dnds <- mouse_dnds %>%
  mutate(omega = dN / dS) %>%
  filter(!is.na(omega) & omega < 10) 

dnds_mouse_gene <- mouse_dnds %>%
  group_by(Ensembl.Gene.ID) %>%
  summarize(mean_dN = mean(dN, na.rm=TRUE),
            mean_dS = mean(dS, na.rm=TRUE),
            omega = mean(omega, na.rm=TRUE))

#### 2) filter data for bhlh factors only, and then combine into one dataset with filtered dndsvalue

bhlh_genes <- bhlh %>%
  filter(DBD=="bHLH") %>%
  select(Ensembl.ID)

bzip_genes <- bhlh %>%
  filter(DBD=="bZIP") %>%
  select(Ensembl.ID)

bhlh_mouse_dnds <- dnds_mouse_gene %>%
  inner_join(bhlh_genes, by=c("Ensembl.Gene.ID"="Ensembl.ID"))

bzip_mouse_dnds <- dnds_mouse_gene %>%
  inner_join(bzip_genes, by=c("Ensembl.Gene.ID"="Ensembl.ID"))

neurod2mouse <- bhlh_mouse_dnds %>%
  filter(Ensembl.Gene.ID == "ENSG00000171532")
neurod2mouse ### dn/ds value of NeuroD2 from mouse versus human



#### 3) Plot

box_bhlh <- bhlh_mouse_dnds %>%
  select(Ensembl.Gene.ID, omega) %>%
  mutate(set = "bHLH TFs")

box_genome <- genome %>%
  select(Ensembl.Gene.ID, omega) %>%
  mutate(set = "Genome-wide")

box_bzip <- bzip_mouse_dnds %>%
  select(Ensembl.Gene.ID, omega) %>%
  mutate(set = "bzip TFs")

box_df <- bind_rows(box_bhlh, box_genome, box_bzip)
neuro_mouse <- as.numeric(neurod2mouse$omega[1])

##side by side violin plot
neuro <- as.numeric(neurod2mouse$omega[1])

ggplot(box_df, aes(x = set, y = omega, fill = set)) +
  geom_violin(alpha = 0.7, trim = TRUE, color = NA) +
  geom_boxplot(width = 0.10, outlier.shape = NA) +
  geom_point(aes(x = set, y = neuro), color = "red", size = 3) +
  scale_fill_manual(values = c(
    "Genome-wide" = "#f79f79",
    "bHLH TFs"    = "#F7D08A",
    "bZIP TFs"    = "#DFF081"
  )) +
  scale_color_manual(values = c(                       # box outlines match fills
    "Genome-wide" = "#f79f79",
    "bHLH TFs"    = "#F7D08A",
    "bZIP TFs"    = "#DFF081"
  )) +
  labs(title = "dN/dS (W) — bHLH & bZIP TFs vs Mouse Genome",
       x = NULL, y = "W (dN/dS)") +
  coord_cartesian(ylim = c(0, 3)) +                 # <- y max = 5
  theme_minimal(base_size = 12) +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))
