#!/bin/bash

cd data

[ -d rrna_16S_ffns ] || mkdir rrna_16S_ffns

for p_genome in genomes/*.fna.gz ; do

  re="(GCA_[0-9]+\.[0-9]+).*$"
  [[ $p_genome =~ $re ]] && genome=${BASH_REMATCH[1]}
  echo $genome
  [[ -e rrna_16S_beds/${genome}_rrna_16S.bed ]] || continue
  gunzip $p_genome
  bedtools getfasta -s -name \
    -fi ${p_genome%.gz} \
    -bed rrna_16S_beds/${genome}_rrna_16S.bed \
    -fo rrna_16S_ffns/${genome}_rrna_16S.ffn
  gzip ${p_genome%.gz}

done

cat rrna_16S_ffns/*.ffn > all_rrna_16S.ffn
