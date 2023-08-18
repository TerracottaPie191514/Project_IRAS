#!/bin/bash
set -e
# concatenates the separate lanes of illumina reads per sample into one sample file for both reads
indir="$HOME/results/00_raw_reads"
outdir="$HOME/results/01_merged_reads/"
mkdir -p $outdir

for dir in $indir/*/
  do
    if (($(ls "$dir"*.fq.gz | wc -l) > 2)); then
          date
          echo "More than 2 fq.gz files present in directory"
          ls "$dir" | tee file_path.tab
          grep -Pro '[F,f]irm\_[0-9]{1,2}\_[0-9]{1,2}\_F[A-Z]{3}19070\d{4}\-\da\_[A-Z0-9]{9}' file_path.tab | sort -u | tee files.tab # this regex captures the names of all samples up until the lanes
          sample_file=$(head -n1 ./files.tab)
          sample_name=$(echo "$sample_file" | cut -f 1-4 -d "_");
          sample_rest=$(echo "$sample_file" | cut -f 5 -d "_");
          echo "Concatenating first pair"
          cat ${dir}${sample_file}_*_1.fq.gz > "$outdir"${sample_name}_${sample_rest}_concat_1.fq.gz;
          echo "Concatenating second pair"
          cat ${dir}${sample_file}_*_2.fq.gz > "$outdir"${sample_name}_${sample_rest}_concat_2.fq.gz;
          rm file_path.tab files.tab
    fi
  done
