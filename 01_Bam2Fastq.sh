#!/bin/bash
###Created by RMeller 2022 OCT
## Convert unaligned Bam to Fastq in specific sample directory

# script called using cat samples.txt | parallel -j 8 Script/01_BamtoFastq.sh

# all calls from $PROJECT folder. 

# create directories

PROJECT=""path_to_project_folder"/ANALYSIS"
mkdir $PROJECT/${1}
mkdir $PROJECT/${1}/01_Fastq
indir="$PROJECT/${1}/01_Fastq"
#cd ${indir}

## Move Bam file to in_directory

echo moving ${1} files to folder ${indir}
cp RawBam/${1}_*.bam ${indir}
cd ${indir}
pwd


## now create a loop to read in the bam files and then use picard tools

for i in *.bam
#	$(ls *.bam | rev | cut -c 5- | rev )
           do echo $i
           File=$(basename $i .bam)
           echo $File
           PicardCommandLine SamToFastq I=${i} FASTQ=$File.fq
           gzip $File.fq ${File}.fq.gz
           fastqc ${File}.fq.gz

done
