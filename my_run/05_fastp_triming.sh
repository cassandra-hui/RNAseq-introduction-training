#!/usr/bin/bash
# trimming the mouseMT fastq files

## This doesn't trim anything


## creating output folder, in case it does not exists
mkdir -p 030_p_trim

INPUT_FOLDER=fastq_files

## each job grab a specific line from sampleNames.txt
for SAMPLE in `cat sampleNames.txt`
do
 fastp -i $INPUT_FOLDER/${SAMPLE}.fastq \
              -o 030_p_trim/${SAMPLE}.trimmed.fastq \ 
              --cut_front --cut_tail \
              --cut_window_size 12 \
              --cut_mean_quality 25 \
              2> 030_p_trim/030_l_trim_out.${SAMPLE}.log
done