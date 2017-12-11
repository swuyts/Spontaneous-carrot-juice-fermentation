#!/bin/bash

threads=8

# decent alignment with the ginsi algorithm
mafft --globalpair --maxiterate 1000 --thread $threads data/all_rrna_16S.ffn > data/all_rrna_16S_ginsi.fasta
