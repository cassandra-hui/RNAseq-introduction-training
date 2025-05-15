#!/usr/bin/bash
# script to perform the QC of the mouseMT fastq files

# creating the output folder
mkdir -p 010_d_fastqc/

fastqc -o 010_d_fastqc fastq_files/*.fastq