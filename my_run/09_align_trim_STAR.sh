#!/usr/bin/bash
# mapping trimmed mouseMT reads



mkdir -p 044_d_STAR_map_trimmed

for SAMPLE in `cat sampleNames.txt`
do

 FASTQ_NAME=030_d_trim/${SAMPLE}.trimmed.fastq

 STAR --runThreadN 4 --genomeDir 041_d_STAR_mouseMT_reference \
                  --outSAMtype BAM SortedByCoordinate \
                  --outFileNamePrefix  044_d_STAR_map_trimmed/${SAMPLE}_trimmed. \
                  --quantMode GeneCounts \
                  --readFilesIn $FASTQ_NAME
done