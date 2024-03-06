# genome_visual
## Using gggenomes package to draw images

3 files are required.

```
library(tidyverse)
library(gggenomes)
  
  
main_genes <- read.csv("genes.bed", sep = "\t", header = F)
main_genes$V5 <- NULL
colnames(main_genes) <- c("seq_id","start","end","name","strand")
main_genes$type <- "CDS"
  
# Order of the sequences 
main_seqs <- read.csv("main_seq.csv", sep = ";", header = T)

main_track <- read.csv("track_link.csv", sep = ";")

g <- gggenomes(genes = main_genes, seqs=main_seqs, links = main_track) +
  # Create lines
  geom_seq() +
  # Add genes, with legends
  #geom_gene(aes(fill=name), size = 5) +
  geom_gene(size = 6) +
  geom_gene_tag(aes(label=name), check_overlap = TRUE, vjust = -1) +
  # Add lable to line (left)
  geom_bin_label() +
  geom_link(aes(fill=perc), offset = 0.15) +
  #geom_link(offset = 0.15) +
  scale_fill_viridis_b() + scale_colour_viridis_b()
g
  
ggsave(g, file="Fig1.eps", device="eps")
ggsave(g, file="Fig1.png", device="png")
```
