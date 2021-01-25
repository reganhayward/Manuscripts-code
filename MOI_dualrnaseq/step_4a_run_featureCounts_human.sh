#!/bin/bash

#Bash script to process bam files for FeatureCounts
#Need to adjust for where bam files are finally stored

FILES=final_bams_human/*.bam
COUNTER=0

#Create list of files
for f in $FILES
do
	size=${#f} #get the length
	length=${f:21:$size} #extract just the file name, removing the folder location
	just_name=${length:0:${#length} - 3} #Removing the .bam
	array_names[$COUNTER]=$just_name #add filename to array
	let COUNTER=COUNTER+1
done

#Remove duplicates
sorted_unique_ids=($(echo "${array_names[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

#Store the locations and make necessary folders
mkdir feature_counts_human
BAM_FILES=final_bams_human/
FINAL_FILE=feature_counts_human/
END="bam"
TXT=".txt"
COUNT=".count"

#Loop through unique names
for i in "${sorted_unique_ids[@]}"
do
   bam_files="${BAM_FILES}${i}${END}"
   count_files="${FINAL_FILE}${i}${TXT}"
   count_files_for_R="${FINAL_FILE}${i}${COUNT}"  

	featureCounts -T 28 -a Homo_sapiens.GRCh38.87.gtf -p -C -Q 10 -o $count_files $bam_files
	tail -n +3 $count_files | cut -f1,7 > $count_file_for_R
	echo "${i}: Has been completed"
done