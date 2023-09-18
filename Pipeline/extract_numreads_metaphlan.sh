#!/bin/bash
set -e

indir="$HOME/results/04_mapped_bac_reads_metaphlan"
outfile="$HOME/results/total_mapped_reads_metaphlan.tab"

# This script will extract the total number of mapped reads from MetaPhlan output (with the output name set as profiled_metagenome.txt)

samples=$(ls $indir | grep -Po '[F,f]irm\_[0-9]{1,2}\_[0-9]{1,2}\_F[A-Z]{3}19070\d{4}\-\da\_[A-Z0-9]{9}' | cut -f 1-3 -d "_")
echo $samples

for sample in $samples
  do
    echo -en $sample'\t' >> $outfile # appends sample without newline character
    sed '7q;d' $indir/*$sample*/profiled_metagenome.txt | cut -f 5 >> $outfile # appends number of mapped reads
  done