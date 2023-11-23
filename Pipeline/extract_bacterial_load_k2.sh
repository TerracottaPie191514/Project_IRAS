#!/bin/bash
set -e

indir="$HOME/results/05_mapped_bac_reads_kraken2/total"
outfile="$HOME/results/bacterial_load_kraken2.tab"

# This script will extract the percentage and total number of bacterial mapped reads from Kraken2 output (with the output name set as *_report.txt)

echo -e "Sample_Unique	ReadPerc	ReadTot" > $outfile #add header line

samples=$(ls $indir | grep -Po '[F,f]irm\_[0-9]{1,2}\_[0-9]{1,2}\.kreport' | cut -f 1 -d ".")
#echo $samples

for sample in $samples
  do
    echo -en $sample'\t' >> $outfile # appends sample without newline character
    sed '4q;d' $indir/*$sample.kreport | cut -f 1,2 >> $outfile # appends number of mapped reads
  done