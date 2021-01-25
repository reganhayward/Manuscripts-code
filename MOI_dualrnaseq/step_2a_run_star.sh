#!/bin/bash

#Bash script to align reads to Host genome using STAR


#Select all files that are .fq.gz in the folder from fastqc
FILES=/trim_galore_files/*.fastq.gz

COUNTER=0
#Create list of files (excluding folder location)
for f in $FILES
do
	size=${#f}
	length=${f:19:$size} #extract just the file name, removing the folder location
	just_name=${length:0:${#length} - 9} #Removing the _r1.fq.gz
	array_names[$COUNTER]=$just_name #add filename to array
	let COUNTER=COUNTER+1
done

#Remove duplicates and sort filenames
sorted_unique_ids=($(echo "${array_names[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

#Store the locations and make necessary folders
RAW_FILES=/trim_galore_files/
ALIGNED_FILES=/star_aligned/
R1="R1.fastq.gz"
R2="R2.fastq.gz"

#Create index
STAR --runThreadN 35 --runMode genomeGenerate --genomeDir STAR_Genome --genomeFastaFiles Homo_sapiens.GRCh38.fa --sjdbGTFfile Homo_sapiens.GRCh38.87.gtf

#Loop through unique names
for i in "${sorted_unique_ids[@]}"
do
   pair1="${RAW_FILES}${i}${R1}"
   pair2="${RAW_FILES}${i}${R2}"
   filename="${ALIGNED_FILES}${i}.bam"
   #Align
   STAR --runThreadN 25 --runMode alignReads --readFilesIn $pair1 $pair2 --readFilesCommand zcat --genomeDir /STAR_Genome --outFileNamePrefix $filename --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 25 --outTmpDir temp/ 
   echo "${i}: Has been completed"
done