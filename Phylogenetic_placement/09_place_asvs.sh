#!/bin/bash

threads=8

# raxml 8.2.11 (June 2017)
raxml=/media/harddrive/tools/standard-RAxML/raxmlHPC-PTHREADS-AVX

cd data

[ -d tree_placement ] || mkdir tree_placement
cd tree_placement
rm *.placement

# convert tsv file with ASV sequences to fasta
tail -n +2 ../../input/ASVs_CJ.tsv | awk -v FS='\t' -v OFS='' '{ 
  print ">", gensub(" ", "_", "g", $4)
  print $1
}' > asvs.fasta

# align ASVs to 16S reference alignment
mafft --addfragments asvs.fasta \
  --reorder --thread $threads \
  ../all_rrna_16S_ginsi.fasta > alignment_with_asvs.fasta

# perform placement of all ASVs
$raxml -T $threads -f v \
  -R RAxML_binaryModelParameters.PARAMS \
  -m GTRGAMMA \
  -s alignment_with_asvs.fasta \
  -t RAxML_result.PARAMS \
  -n placement \
  --epa-keep-placements=100 \
  --epa-prob-threshold=0.5

