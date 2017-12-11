#!/bin/bash

# raxml 8.2.11 (June 2017)
raxml=/media/harddrive/tools/standard-RAxML/raxmlHPC-PTHREADS-AVX

cd data

[ -d tree_16S ] || mkdir tree_16S
cd tree_16S

$raxml -T 8 -f a \
  -m GTRCAT \
  -p 1991 \
  -x 1991 -N autoMRE \
  -s ../all_rrna_16S_ginsi.fasta \
  -n lgc_16S
