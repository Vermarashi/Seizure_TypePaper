#!/bin/bash
#samtools merge and manipulate
#adapted from RMeller 1/12/2022

echo "running on $HOSTNAME"
echo "sample name = $1"
start=$(date +%s)
echo "start time: $start"

## Set up directories
PROJECT="path_to_project_folder"

# make a results folder and designate as output
mkdir ${PROJECT}/ANALYSIS/${1}/05_Results
out_dir=${PROJECT}/ANALYSIS/$1/05_Results

# Now move to the aligned data and set as indirectory
in_dir2=${PROJECT}/ANALYSIS/$1/03_GRCH38_AlignedData
cd ${in_dir2}

echo "samtools version: "
samtools --version


#Calls
# create file with name of all of the bam files in it
call1a="ls */*Aligned.out.bam > bam.txt"
echo $call1a
eval $call1a

# Merge files
call1="samtools merge -@ 12 -r -b bam.txt $out_dir/$1_Combined_Aligned.out.bam"
# Sort Files
call2="samtools sort -m 2500M --threads 12 -o $out_dir/$1_Combined_Aligned.SortedByCord.out.bam $out_dir/$1_Combined_Aligned.out.bam"
# index files
call3="samtools index -@ 12 $out_dir/$1_Combined_Aligned.SortedByCord.out.bam "

echo $call1
eval $call1

echo $call2
eval $call2

echo $call3
eval $call3

end=$(date +%s)
echo "end time: $end"
runtime_s=$(echo $(( end - start )))
echo "total run time(s): $runtime_s"
sec_per_min=60
sec_per_hr=3600
runtime_m=$(echo "scale=2; $runtime_s / $sec_per_min;" | bc)
echo "total run time(m): $runtime_m"
runtime_h=$(echo "scale=2; $runtime_s / $sec_per_hr;" | bc)
echo "total run time(h): $runtime_h"
