#!/usr/bin/Rscript

library(tidyverse)
library(stringr)

if (! dir.exists("data/ftp_paths")) dir.create("data/ftp_paths")

genomes_genbank <- "ftp://ftp.ncbi.nih.gov/genomes/genbank/assembly_summary_genbank.txt" %>%
  read_tsv(file = ., skip = 1, na = c("", "na"), quote = "")  %>%
  mutate(wgs_master2 = str_extract(wgs_master, "^[^.]+")) %>%
  rename(assembly_accession = `# assembly_accession`)

genomes_sunetal <- read_csv2(file = "input/sunetal_strainlist_frompdf.csv") %>%
  select(species = `Species Name`, strain = `StrainID`, wgs_master2 = `Genome Accession`) %>%
  left_join(genomes_genbank) %>%
  select(species, strain, wgs_master2, assembly_accession)

# Unreadable genomes from PDF
difficult_genomes <- c(
  "NC_008530" = "GCA_000014425.1",
  "CP003851 - CP00385" = "GCA_000300135.1",
  "Q489736 - DQ48974" = "GCA_000026405.1",
  "FN822744" = "GCA_000196855.1",
  "CP001753 - CP00175" = "GCA_000092505.1",
  "C_008496, NC_0085" = "GCA_000014445.1",
  "C2KK01" = "GCA_000160595.1",
  "JQBE00000000" = "GCA_001640785.1"
)

genomes_sunetal <- genomes_sunetal %>%
  mutate(assembly_accession = ifelse(is.na(assembly_accession), difficult_genomes[wgs_master2], assembly_accession))

genomes_sunetal %>%
  select(assembly_accession, species, strain) %>%
  write_tsv(col_names = T, path = "data/genomes.tsv")

genomes_sunetal %>%
  select(assembly_accession) %>%
  left_join(genomes_genbank) %>%
  mutate(link = str_c(ftp_path, assembly_accession, sep = "/")) %>%
  mutate(link = str_c(link, asm_name, sep = "_")) %>%
  select(link) %>%
  write_tsv(col_names = F, path = "data/ftp_paths/ftp_paths.txt")
