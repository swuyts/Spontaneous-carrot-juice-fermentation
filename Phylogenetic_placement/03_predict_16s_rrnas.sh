#!/bin/bash

cd data

[ -d rrna_gffs ] || mkdir rrna_gffs
[ -d rrna_16S_gffs ] || mkdir rrna_16S_gffs

for p_genome in genomes/*.fna.gz ; do

  gunzip $p_genome
  re="(GCA_[0-9]+\.[0-9]+).*$"
  [[ $p_genome =~ $re ]] && genome=${BASH_REMATCH[1]}
  echo $genome
  ../tools/barrnap/bin/barrnap ${p_genome%.gz} > rrna_gffs/${genome}_rrnas.gff
  grep 'Name=16S_rRNA' rrna_gffs/${genome}_rrnas.gff > rrna_16S_gffs/${genome}_rrna_16S.gff
  gzip ${p_genome%.gz}

done

