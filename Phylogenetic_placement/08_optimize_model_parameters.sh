#!/bin/bash

# raxml 8.2.11 (June 2017)
raxml=/media/harddrive/tools/standard-RAxML/raxmlHPC-PTHREADS-AVX

cd data

[ -d tree_placement ] || mkdir tree_placement
cd tree_placement

$raxml -T 8 -f e \
  -m GTRGAMMA \
  -s ../all_rrna_16S_ginsi.fasta \
  -t ../tree_16S/RAxML_bestTree.lgc_16S \
  -n PARAMS

