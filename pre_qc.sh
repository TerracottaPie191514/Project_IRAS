#!/bin/bash
set -e

indir="$HOME/results/00_raw_reads/"
mergedir="$HOME/results/01_merged_reads/"
outdir="$HOME/results/QC/untrimmed_reads/"
mkdir -p $outdir

echo "QC on merged files"

for file in $mergedir*.gz
  do
    echo $file
    date
    cd $outdir
    # We have to copy the original files due to permission issues
    cp $file .
    echo "Running FastQC"
    fastqc *.gz
    rm *.gz
  done

echo "QC on non-merged files"

for dir in $indir/*/
  do
    if (($(ls "$dir"*.fq.gz | wc -l) < 3)); then
      echo "less than 3 fq.gz files in directory"
      cd $outdir
      cp $dir/*.gz .
      echo "Running FastQC"
      fastqc *.gz
      rm *.gz
    fi
  done

echo "Unzipping..."
for filename in *.zip
  do
  unzip $filename
done

echo "Creating summary"
cat */summary.txt > ./fastqc_summaries.txt
multiqc .
