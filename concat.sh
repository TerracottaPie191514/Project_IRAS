#!/bin/bash
set -e
INDIR="$HOME/results/03_mapped_reads"
OUTDIR="${INDIR}_combined"
OUTDIR2="${OUTDIR}_selectColumns"
FINALOUTFILE="$HOME/results/Resfinder_mapped_AllSamples_Longformat.tab"

mkdir -p $OUTDIR
mkdir -p $OUTDIR2

date

samples=$(ls $INDIR | grep -Po '[F,f]irm\_[0-9]{1,2}\_[0-9]{1,2}\_F[A-Z]{3}19070\d{4}\-\da\_[A-Z0-9]{9}' | cut -f 1-3 -d "_")
echo $samples


for sample in $samples
 do
 	# create new sample all files
 	echo $sample
 	echo -e "name	%unambiguousReads	unambiguousMB	%ambiguousReads	ambiguousMB	unambiguousReads	ambiguousReads	assignedReads	assignedBases" > $OUTDIR/$sample\_all.tab #add header line without #

 	# list all scafstats in various subfolders of sample
 	# filter out all header lines using inverse grep matching -v
 	cat $INDIR/*$sample*/scafstats.txt | grep -v -P "^#" >> $OUTDIR/$sample\_all.tab
 done

echo -e "\n\nProcessing two relevant columns and adding sample name for concat later\n"

for file in $OUTDIR/*.tab
 do
 	sample=`basename $file _all.tab`
 	echo $sample

 	# add sample using bash variable
 	# restore header line which got also actual sample instead "Sample" column name
 	awk -v sample="$sample" '{ print sample "\t" $1 "\t" $8 }' $file \
 			| sed -r -e "s/^\S+\tname\t(.+)/Sample\tARG\t\1/" \
 			> $OUTDIR2/${sample}_columns.tab
 done
 
echo -e "\n\nNow combine all files from all samples in a single large table\n\n"
date

echo -e "SampleName\tARG\tAssignedReads" > $FINALOUTFILE

for file in $OUTDIR2/*.tab
 do
 	echo $file
 	# print content each file except the top header line.
 	tail -n +2 $file >> $FINALOUTFILE
 done
