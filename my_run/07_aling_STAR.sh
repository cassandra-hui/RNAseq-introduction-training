#!/usr/bin/bash
# aligning mouseMT reads with STAR

###############################################################
# Had to mannualy install and older version and use it. There is a bug for Macs
# Go to the fodlder

cd /Users/cassandrahui/Documents/Programs/star-2.7.5a-0

# Export the path
# While in the star-2.7.5a-0 directory
export PATH="$PWD/bin:$PATH"
#################################################################


mkdir -p 042_d_STAR_mouseMT_map_raw

for SAMPLE in `cat sampleNames.txt`
do
 FASTQ_NAME=fastq_files/${SAMPLE}.fastq

 STAR --runThreadN 4 --genomeDir 041_d_STAR_mouseMT_reference \
               --outSAMtype BAM SortedByCoordinate \
               --outFileNamePrefix  042_d_STAR_map_raw/${SAMPLE}. \
               --quantMode GeneCounts \
               --readFilesIn $FASTQ_NAME 
done