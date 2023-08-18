#!/bin/bash
set -e

indir="$HOME/results/01_merged_reads"
indir2="$HOME/results/00_raw_reads"
outdir="$HOME/results/02_cleaned_reads"
mkdir -p $outdir

echo "Clean-up on merged files"

for file in $indir/*_1.fq.gz
 do
    date
    echo $file
    sample=`basename $file _1.fq.gz`
    
    file1=$indir/$sample\_1.fq.gz
    file2=$indir/$sample\_2.fq.gz
    fileout1=$outdir/$sample\_1_clean.fq.gz
    fileout2=$outdir/$sample\_2_clean.fq.gz

    /home/melissen/miniconda3/bin/bbduk.sh   \
        in1=$file1 \
        in2=$file2 \
        out1=$fileout1 \
        out2=$fileout2 \
        ref=/home/melissen/miniconda3/opt/bbmap-39.01-1/resources/adapters.fa \
        k=21 mink=6 ktrim=r ftm=5 qtrim=rl \
        trimq=20 minlen=30 \
        overwrite=true -Xmx800m
        
 done

echo "Clean-up on non-merged files"

for dir in $indir2/*/
  do
    if (($(ls "$dir"*.fq.gz | wc -l) < 3)); then
      echo "less than 3 fq.gz files in directory"
      for file in $dir*_1.fq.gz
        do 
          date
          echo $file
          sample=`basename $file _1.fq.gz`
          cd $outdir
          # We have to copy the original files due to permission issues
          cp $dir/*.fq.gz .
          
          echo "running BBduk"
          file1=$outdir/$sample\_1.fq.gz
          file2=$outdir/$sample\_2.fq.gz
          fileout1=$outdir/$sample\_1_clean.fq.gz
          fileout2=$outdir/$sample\_2_clean.fq.gz
      
          /home/melissen/miniconda3/bin/bbduk.sh   \
              in1=$file1 \
              in2=$file2 \
              out1=$fileout1 \
              out2=$fileout2 \
              ref=/home/melissen/miniconda3/opt/bbmap-39.01-1/resources/adapters.fa \
              k=21 mink=6 ktrim=r ftm=5 qtrim=rl \
              trimq=20 minlen=30 \
              overwrite=true -Xmx800m
            # Remove the original *.fq.gz files but not the *_clean.fq.gz files
            find . -type f -name "*.fq.gz" ! -name "*_clean.fq.gz" -delete
        done
    fi
  done

echo -e "\n\n------ Cleanup done ------\n\n"
date
