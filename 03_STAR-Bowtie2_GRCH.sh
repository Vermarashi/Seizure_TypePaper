#!/bin/bash
###Again run this with a cat samples.txt | parallel ./star.sh
#######################################################################
# STAR  w/STAR-FUSION Script and Bowtie2 for unmatched reads SEIZURES

echo "running on $HOSTNAME"
echo "sample name = $1"

## make and define folders
PROJECTS="path_to_project_folder"

# create output folder
mkdir ${PROJECTS}/ANALYSIS/${1}/03_Grch38_AlignedData
out_dir="${PROJECTS}/ANALYSIS/${1}/03_Grch38_AlignedData"

# now define input folder
in_dir="${PROJECTS}/ANALYSIS/${1}/02_TrimmedFastq"
cd ${in_dir}

## define files for STAR
REF="/REF/REFERENCES/GRCH38_MainChrs"
echo ${REF}

GTF="/REF/REFERENCES/GRCH38_MainChrs/Homo_sapiens.GRCh38.109.gtf"
echo ${GTF}

BOWTIE2REF="/REF/REFERENCES/GRCH38_MainChrs/BOWTIE2_INDEX/GRCH38_MAINChrs"
echo ${BOWTIE2REF}


# You will have a bunch of fastq.gz files.  First extract the core name.
# Remove the extras at the end of the name, so that the resultant bam file looks nice

for fastq in $(ls *_trimmed.fq.gz | rev | cut -c 15-| rev )

do
	echo "STAR version: "
	STAR --version
	mkdir ${out_dir}/${fastq}
	echo "SAMPLE: ${fastq}"

#This creates an item (call of the program) the call is echo’ed and then ran with eval. 
# this is useful for running large numbers of files with and collecting log files for the run, so you can see what went wrong. 
 

	call1="STAR --twopassMode Basic \
	--limitBAMsortRAM 20000000000 \
	--genomeDir ${REF} \
        --outReadsUnmapped Fastx \
	--sjdbGTFfile ${GTF} \
	--outFilterType BySJout \
	--outSAMattributes NH HI AS nM NM MD jM jI MC ch \
	--outSAMattrRGline ID:${fastq} PL:ION-TORRENT PU:${fastq} LB:${1} SM:${1} \
	--outFilterMultimapNmax 5 \
	--outFilterMismatchNmax 999 \
	--outFilterMismatchNoverReadLmax 0.04 \
	--alignIntronMin 10 \
	--alignIntronMax 100000 \
	--alignMatesGapMax 100000 \
	--alignSJoverhangMin 8 \
	--alignSJDBoverhangMin 1 \
	--sjdbScore 1 \
	--readFilesCommand zcat \
	--runThreadN 8 \
	--chimOutType Junctions SeparateSAMold WithinBAM SoftClip \
	--chimOutJunctionFormat 1 \
	--chimSegmentMin 15 \
	--limitOutSJcollapsed 5000000 \
	--outSAMtype BAM SortedByCoordinate \
	--quantMode TranscriptomeSAM GeneCounts \
	--outSAMheaderHD @HD VN:1.4 SO:coordinate \
	--outFileNamePrefix ${out_dir}/${fastq}/${fastq} \
	--readFilesIn ${in_dir}/${fastq}_trimmed.fq.gz"

 	echo ${call1}
	eval ${call1}

	call2="bowtie2 --local --very-sensitive-local -p 8 -q --mm \ 
	-x $BOWTIE2REF -U ${out_dir}/${fastq}/${fastq}Unmapped.out.mate1 \
	 --un ${out_dir}/${fastq}/${fastq}_sbt2_unmap.fq | samtools view -bS - | samtools sort -o ${out_dir}/${fastq}/${fastq}_unmapped_remapped_remapBowtie2.bam - "
	 
	echo ${call2}
	eval ${call2}

	call3="PicardCommandLine MergeSamFiles USE_THREADING=true MSD=true AS=true \
	I=${out_dir}/${fastq}/${fastq}Aligned.sortedByCoord.out.bam \
	I=${out_dir}/${fastq}/${fastq}_unmapped_remapped_remapBowtie2.bam \
	O=${out_dir}/${fastq}/${fastq}_STARBowtie2.bam"
	
	echo ${call3}
	eval ${call3}

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
