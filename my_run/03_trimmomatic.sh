#!/usr/bin/bash
# trimming the mouseMT fastq files


## creating output folder, in case it does not exists
mkdir -p 030_d_trim

INPUT_FOLDER=fastq_files

## each job grab a specific line from sampleNames.txt
for SAMPLE in `cat sampleNames.txt`
do
 trimmomatic SE -phred33 \
              $INPUT_FOLDER/${SAMPLE}.fastq \
              030_d_trim/${SAMPLE}.trimmed.fastq \
              ILLUMINACLIP:/data/DATA/adapters/TruSeq3-PE.fa:2:30:10 \
              SLIDINGWINDOW:3:25 2> 030_d_trim/030_l_trim_out.${SAMPLE}.log
done