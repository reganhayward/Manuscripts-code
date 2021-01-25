#!/bin/bash

#-------------
#Removing reads from other species
#-------------


#The order
#1)  Remove duplicates
#2)  Sort by readname, helps extract the paired end reads
#3)  Extract r1 and r2 using bedtools
#4)  Running Bowtie to align
#5)  Extracting only mapped reads
#6)  Convert only mapped to fq  
#7)  Change headers in fq files
#8)  Removing the other species reads from the original BAM
#9)  Getting mapping stats of final bams



#Capture the output
exec > >(tee "processing_out_human_reads_from_chlamydial_samples.log") 2>&1


#---
#1)  Remove duplicates
#---
echo "#---------------"
echo "#1/9) Removing duplicates"
echo "#---------------"

mkdir removed_duplicates_chlamydial

#find/store bam files
bam_files="find *.bam -maxdepth 1 ! -type d"
#loop through
for i in $bam_files
do
   filename_no_ext==${i::-4}
   filename="${filename_no_ext}_no_dups.bam"
   samtools view -b -f 2 -F 1024 -@ 28 ${i} > removed_duplicates_chlamydial/${filename}
   echo "${i}: Has been completed"
done

echo " "
echo "----------------------"
echo " "
echo "#Removing duplicates - done"
echo " "




#---
#2)  Sort by readname, helps extract the paired end reads
#---
echo "#---------------"
echo "#2/9) Sorting BAMs by name"
echo "#---------------"


mkdir bam_sorted_name_chlamydial

#find/store bam files
bam_files="find removed_duplicates_chlamydial/*.bam -maxdepth 1 ! -type d"
#loop through
for i in $bam_files
do
   filename_no_ext==${i::-4}
   filename="${filename_no_ext}_sorted_name.bam"
   samtools sort -@ 28 -n removed_duplicates/${i} > bam_sorted_name_chlamydial/${filename}
   echo "${i}: Has been completed"
done

echo " "
echo "----------------------"
echo " "
echo "#Sorting BAMs by name - done"
echo " "

#---
# Steps 3:9)  Extract r1 and r2 using bedtools, then align to host genome using Bowtie2
#---

echo "#---------------"
echo "# Running Steps 3:9 "
echo "#---------------"

mkdir r1_r2
mkdir bowtie_mapped_human
mkdir bowtie_mapped_human_only_mapped
mkdir final_bams
mkdir final_bams_stats

bowtie2-build Homo_sapiens.GRCh38.fa Human_genome/Homo_sapiens.GRCh38
GENOME_INDEX="Human_genome/Homo_sapiens.GRCh38"
THREADS=28

#find/store bam files
bam_files="find bam_sorted_name_chlamydial/*.bam -maxdepth 1 ! -type d"
#loop through
for i in $bam_files
do
   filename_simple=${i::-24}
   filename_simple_r1="${filename_simple}_r1.fq"
   filename_simple_r2="${filename_simple}_r2.fq"
   filename="${filename_no_ext}_sorted_name.bam"
   filename_ct="${filename_simple}_human_mapped.bam"
   filename_ct_fq="${filename_simple}_human_mapped.fq"
   filename_final="${filename_simple}_human_removed.bam"
   filename_stats="${filename_simple}_bamstats.log"
   
   #extract .fq (step 3)
   bedtools bamtofastq -i bam_sorted_name_chlamydial/${i} -fq r1_r2/${filename_simple_r1} -fq2 r1_r2/${filename_simple_r2}
   #align (step 4)
   bowtie2 -p $THREADS -x $GENOME_INDEX --phred33 -N 1 --very-sensitive-local -1 r1_r2/${filename_simple_r1} -2 r1_r2/${filename_simple_r2} \
    | samtools view -@ $THREADS -bS -h -F4 - > bowtie_mapped_human/${filename_ct}
   #Extract only mapped reads (to Human) (step 5)
   samtools view -@ $THREADS -bSh -F4 bowtie_mapped_human/${filename_ct} -o bowtie_mapped_human_only_mapped/${filename_ct}
   #Convert all the human mapped reads into fq file (step 6)
   bedtools bamtofastq -i bowtie_mapped_human_only_mapped/${filename_ct} -fq bowtie_mapped_human_only_mapped/${filename_ct_fq}
   #For some reason, an @ gets put before the SN, so have to remove it (step 7)
   sed -i -e 's/@SN/SN/g' bowtie_mapped_human_only_mapped/${filename_ct_fq}
   #Removing the other species reads from the original BAM" (step 8)
   picard FilterSamReads I=${i} O=final_bams/${filename_final} READ_LIST_FILE=bowtie_mapped_human_only_mapped/${filename_ct_fq} FILTER=excludeReadList
   #Saving down the final bam stats" (step 9)
   bamtools stats -in final_bams/${filename_final} > final_bams_stats/${filename_stats}

   echo "${i}: Has been completed"

done


echo " "
echo "----------------------"
echo " "
echo " Finished"
echo " "