#!/bin/bash

# Trim Galore Script
# to be used with cat samples.txt| parallel -j n /path/to/./trimgalore.sh

# create directories
PROJECT="/Path_to _project_folder/"


indir="$PROJECT/ANALYSIS/${1}/01_Fastq"
mkdir $PROJECT/ANALYSIS/${1}/02_TrimmedFastq
outdir="$PROJECT/ANALYSIS/${1}/02_TrimmedFastq"
cd $indir
pwd

## Now loop trimgalore on *.fq.gz files in directory
for fastq in  *.fq.gz   #use for pairs#$( ls *_.fq | rev | cut -c -6 | rev | uniq)
	do
	echo ${fastq}
	trim_galore --fastqc -a ATCACCGACTGCCCATAGAGAGGCTGAGAC -q 20 --gzip ${fastq} --output_dir ${outdir}
done


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
