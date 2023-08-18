#!/bin/bash
set -e

indir="$HOME/results/02_cleaned_reads"
outdir="$HOME/results/QC/trimmed_reads_1"
mkdir -p $outdir

echo "QC on trimmed reads"

for file in $indir/*.gz
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


echo "Unzipping..."
for filename in *.zip
  do
  unzip $filename
done

echo "Creating summary"
cat */summary.txt > ./fastqc_summaries.txt
multiqc .
