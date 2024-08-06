#!/bin/bash

####Rashi Verma 26 July 2024, Version = Stringtie2countMatrix_v01.sh
####Script to genetrate count matrix for genes and transcripts from bam files using stringtie
####Script will use bam_files as input and run it through stringtie which assembled_gtf -> merged_gtf -> estimate_abundance.
####Finaly produces Gene_count_matrix and Transcript_count_matrix through prepDE.py script

start=$(date +%s)
echo "start time: $start"
echo "running on $HOSTNAME"


# Define the path to your BAM files directory
BAM_DIR="bam_files/"
# Define the path to your GTF file
GTF_FILE="Homo_sapiens.GRCh38.110.gtf"
# Define the output directory
OUTPUT_DIR="analysis/"
# Create the output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"


# Initialize an array to store the GTF file paths
GTF_FILES=()
# Define a function to process each BAM file
assemble_gtf() {
    local BAM_FILE="$1"
    local OUTPUT_DIR="$2"
    local GTF_FILE="$3"
    
    # Get the filename without extension (prefix)
    local PREFIX=$(basename "${BAM_FILE}" .bam)

    # Extract the group name (e.g., G1, G2, G3) from the prefix
    GROUP="${PREFIX%%_*}"

    # Define the output GTF filename with the group suffix
    local OUTPUT_GTF="${OUTPUT_DIR}${PREFIX}.gtf"
    
    # Run stringtie with the specified parameters
    stringtie "${BAM_FILE}" -l "${PREFIX}" -p 8 -G "${GTF_FILE}" -o "${OUTPUT_GTF}"
    
    echo "Generated ${OUTPUT_GTF} for group ${GROUP}"
}
# Export the function so it's available to parallel
export -f assemble_gtf
# Iterate over each BAM file in the directory and run in parallel
find "${BAM_DIR}" -type f -name "*.bam" | parallel -j8 assemble_gtf {} "${OUTPUT_DIR}" "${GTF_FILE}"
# Create the mergelist.txt file by listing all GTF files in the output directory
ls "${OUTPUT_DIR}"*.gtf > "${OUTPUT_DIR}mergelist.txt"
echo "Created mergelist.txt in ${OUTPUT_DIR}"


# Run stringtie --merge using the mergelist.txt file
stringtie --merge -p 8 -G "${GTF_FILE}" -o "${OUTPUT_DIR}stringtie_merged.gtf" "${OUTPUT_DIR}mergelist.txt"
echo "Generated ${OUTPUT_DIR}stringtie_merged.gtf"


# Define the abundance directory
ABUNDANCE_DIR="abundance/"
# Create the abundance directory if it doesn't exist
mkdir -p "${ABUNDANCE_DIR}"
# Define a function to run stringtie for abundance estimation
estimate_abundance() {
    local BAM_FILE="$1"
    local OUTPUT_DIR="$2"
    local ABUNDANCE_DIR="$3"
    
    # Get the filename without extension (prefix)
    local PREFIX=$(basename "${BAM_FILE}" .bam)
    
    # Define the subfolder path within the abundance directory
    local SUBFOLDER="${ABUNDANCE_DIR}${PREFIX}/"
    
    # Create the subfolder
    mkdir -p "${SUBFOLDER}"
    
    # Run stringtie with -e -B and output to the subfolder
    stringtie -e -B -p 8 -G "${OUTPUT_DIR}stringtie_merged.gtf" -o "${SUBFOLDER}${PREFIX}.gtf" "${BAM_FILE}"
    
    echo "Generated ${SUBFOLDER}${PREFIX}.gtf"
}
# Export the function for use with parallel
export -f estimate_abundance
# Use parallel to run stringtie for abundance estimation in parallel
find "${BAM_DIR}" -type f -name "*.bam" | parallel -j8 estimate_abundance {} "${OUTPUT_DIR}" "${ABUNDANCE_DIR}"


# Create and populate count.txt by listing GTF files in the abundance directory
COUNT_FILE="count.txt"
# List GTF files in the abundance directory and add them to count.txt
find "${ABUNDANCE_DIR}" -type f -name "*.gtf" | while read -r GTF_PATH; do
    # Extract the prefix from the path
    PREFIX=$(basename "$(dirname "${GTF_PATH}")")
    
    # Append to count.txt using the desired path format
    echo -e "${PREFIX}\t${GTF_PATH}" >> "${COUNT_FILE}"
done
echo "Created ${COUNT_FILE}"


# Run ./prepDE.py script with count.txt as input
./prepDE.py -i "${COUNT_FILE}" -g gene_count_matrix.csv -t transcript_count_matrix.csv
echo "Finished running ./prepDE.py"


end=$(date +%s)
echo "end time: $end"
runtime_s=$(echo $(( end - start )))
echo " time in sec: $runtime_s"
sec_per_min=60
sec_per_hr=3600
runtime_m=$(echo "scale=2; $runtime_s / $sec_per_min;" | bc)
echo "total run time Min: $runtime_m"
runtime_h=$(echo "scale=2; $runtime_s / $sec_per_hr;" | bc)
echo "total run time hours : $runtime_h"
