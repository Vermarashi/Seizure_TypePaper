# Seizure_2
##Differential Transcript Expression and RNA Splicing to Distinguish Seizure Type
These scripts are to be used on data present in the paper. 

To successfully run this script, you will need the following input files organized in one folder:
 
1.	Bam_files: This folder should contain the BAM files that are essential for the analysis.
2.	Python_script (prepDE.py) (download from https://github.com/gpertea/stringtie/blob/master/prepDE.py)
3.	Stringtie2countMatrix.sh: This shell script is integral to the initial steps of the analysis and will generate the gene_count_matrix and transcript_count_matrix (attached).
 
To get started, you can follow these steps:
Begin by running the "Stringtie2countMatrix.sh" script. This will generate the necessary count matrices using prepDE.py script, including the gene_count_matrix and transcript_count_matrix.
