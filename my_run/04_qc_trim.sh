#!/usr/bin/bash
# QC abalysis of the mouseMT trimmed reads

## fastQC on trimmed fastq files
fastqc 030_d_trim/*.fastq -o 030_d_trim


## multiqc on the fastQC reports AND the trimmomatic logs
multiqc -n 032_r_multiqc_mouseMT_trimmed.html -f --title trimmed_fastq 030_d_trim/
