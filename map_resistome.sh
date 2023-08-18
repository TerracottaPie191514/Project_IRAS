#!/bin/bash
set -e

indir="$HOME/results/02_cleaned_reads"
outdir="$HOME/results/03_mapped_reads"
mkdir -p $outdir


resfinderdb="$HOME/Resfinder/Resfinder_20200127_all_genes.uniq.fasta"


for file in $indir/*_1_clean.fq.gz
 do
 	date
 	echo $file
 	sample=`basename $file _1_clean.fq.gz`
 	infile1=$file
 	infile2=$indir/$sample\_2_clean.fq.gz
  echo $indir/$infile1
 	OUTDIR=$outdir/$sample
	# do actual mapping in SEMIperfectModus
	# disabled outu=$OUTDIR/unmapped.fq
	/home/melissen/miniconda3/bin/bbmap.sh -Xmx28g \
			 in1=$infile1 in2=$infile2 out=/dev/null ref=$resfinderdb \
			 semiperfectmode=t \
			 outm=/dev/null \
			 covstats=$OUTDIR/constats.txt covhist=$OUTDIR/covhist.txt basecov=$OUTDIR/basecov.txt bincov=$OUTDIR/bincov.txt \
			 bhist=$OUTDIR/bhist.txt qhist=$OUTDIR/qhist.txt aqhist=$OUTDIR/aqhist.txt lhist=$OUTDIR/lhist.txt ihist=$OUTDIR/ihist.txt \
			 ehist=$OUTDIR/ehist.txt qahist=$OUTDIR/qahist.txt indelhist=$OUTDIR/indelhist.txt mhist=$OUTDIR/mhist.txt \
			 gchist=$OUTDIR/gchist.txt idhist=$OUTDIR/idhist.txt scafstats=$OUTDIR/scafstats.txt

 done
