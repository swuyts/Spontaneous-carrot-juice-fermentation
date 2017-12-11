#!/bin/bash

cd data

[ -d genomes ] || mkdir genomes
cd genomes

for ftp_path in $(cat ../ftp_paths/ftp_paths.txt) ; do
  # wget --reject '_cds_' ${ftp_path}/*_genomic.fna.gz
  wget ${ftp_path}_genomic.fna.gz
done
