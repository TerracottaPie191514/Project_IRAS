#!/bin/bash
set -e

indir="$HOME/results/02_cleaned_reads"
outdir="$HOME/results/04_mapped_bac_reads"
mkdir -p $outdir


for file in $indir/*_1_clean.fq.gz
 do
 	date
 	sample=`basename $file _1_clean.fq.gz`
 	echo $sample
 	infile1=$file
 	infile2=$indir/$sample\_2_clean.fq.gz
   	OUTDIR=$outdir/$sample
 	mkdir -p $OUTDIR
  metaphlan $infile1,$infile2 --bowtie2out $OUTDIR/metagenome.bowtie2.bz2 \
  -t rel_ab_w_read_stats --nproc 28 --input_type fastq -o $OUTDIR/profiled_metagenome.txt
 done
