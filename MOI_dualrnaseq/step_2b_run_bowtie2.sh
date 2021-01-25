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
ALIGNED_FILES=/bowtie2_aligned/
R1="R1.fastq.gz"
R2="R2.fastq.gz"
THREADS=28

#Create index
bowtie2-build charm.fa Charm_genome/charm_genome
GENOME_INDEX="Charm_genome/charm_genome"

#Loop through unique names
for i in "${sorted_unique_ids[@]}"
do
   pair1="${RAW_FILES}${i}${R1}"
   pair2="${RAW_FILES}${i}${R2}"
   filename="${ALIGNED_FILES}${i}.sam"
   filename_bam="${ALIGNED_FILES}${i}_sorted_coord.bam"
   #Align
   bowtie2 -p $THREADS -x $GENOME_INDEX --phred33 -N 1 --very-sensitive-local -1 $pair1 -2 $pair2 -S $filename
   #sort and convert to BAM
   samtools view -@ 28 -bS $filename | samtools sort -@ 28 - -o $filename_bam
   
   echo "${i}: Has been completed"
done






