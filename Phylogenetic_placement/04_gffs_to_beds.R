#!/usr/bin/Rscript

library(tidyverse)
library(stringr)

if (! dir.exists("data/rrna_16S_beds")) dir.create("data/rrna_16S_beds")

gff_files <- list.files("data/rrna_16S_gffs/", full.names = T)

for (gff_file in gff_files) {
  
  genome <- str_match(gff_file, "([A-Z]+_[0-9]+\\.[0-9]).*\\.gff")[1, 2]
  
  print(genome)
  
  if (file.info(gff_file)$size == 0) next
  
  read_tsv(gff_file, col_names = F) %>%
    select(chrom = X1, start = X4, end = X5, strand = X7) %>%
    mutate(genome = !! genome) %>%
    mutate(name = str_c(genome, 1:n(), sep = "_")) %>%
    mutate(start = start - 1) %>%
    mutate(score = ".") %>%
    select(chrom, start, end, name, score, strand) %>%
    write_tsv(path = str_c("data/rrna_16S_beds/", genome, "_rrna_16S.bed"), col_names = F)
  
}


