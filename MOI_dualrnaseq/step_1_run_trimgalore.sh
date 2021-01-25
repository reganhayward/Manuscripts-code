#!/bin/bash

#Bash script to process different fastq files for TrimGalore


#Select all files that are .fq.gz in the folder from fastqc
FILES=/MOI_files/*.fastq.gz

COUNTER=0
#Create list of files (excluding folder location)
for f in $FILES
do
	size=${#f}
	length=${f:11:$size} #extract just the file name, removing the folder location
	just_name=${length:0:${#length} - 9} #Removing the _r1.fq.gz
	array_names[$COUNTER]=$just_name #add filename to array
	let COUNTER=COUNTER+1
done

#Remove duplicates and sort filenames
sorted_unique_ids=($(echo "${array_names[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

#Store the locations and make necessary folders
RAW_FILES=/MOI_files/
TRIMMED_FILES=/trim_galore_files/
R1="R1.fastq.gz"
R2="R2.fastq.gz"

#Loop through unique names
for i in "${sorted_unique_ids[@]}"
do
   pair1="${RAW_FILES}${i}${R1}"
   pair2="${RAW_FILES}${i}${R2}"
   filename="${TRIMMED_FILES}${i}_trimgalore.log"
   #Run trim_galore command
   trim_galore -q 14 --paired --length 20 --phred33 --gzip -o trim_galore_files/ $pair1 $pair2 &>$filename
   echo "${i}: Has been completed"
done