#!/usr/bin/bash
# multiQC to group QC reports of the mouseMT fastq files

multiqc -n 020_r_multiqc_mouseMT.html -f --title raw_fastq 010_d_fastqc/