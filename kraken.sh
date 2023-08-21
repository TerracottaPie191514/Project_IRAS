#!/bin/bash
set -e

DB=~/results/k2_standard_20230605
indir="$HOME/results/02_cleaned_reads"
outdir="$HOME/results/05_mapped_bac_reads_kraken2"
mkdir -p $outdir

for file in $indir/*_1_clean.fq.gz
 do
  date
  sample=`basename $file _1_clean.fq.gz`
 	echo $sample
 	infile1=$file
 	infile2=$indir/$sample\_2_clean.fq.gz
 	OUTDIR=$outdir/$sample
  kraken2 -db $DB --threads 24 --report $outdir/"$sample"_report.txt --gzip-compressed \
  --confidence 0.2 --memory-mapping --output - --paired $infile1 $infile2
 done
