#!/usr/bin/bash
# indexing the mouse mitochondrial genome

G_FASTA=ref/Mus_musculus.GRCm39.dna.chromosome.MT.fa
G_GTF=ref/Mus_musculus.GRCm39.114.chr.gtf


mkdir -p 041_d_STAR_mouseMT_reference

STAR --runMode genomeGenerate \
     --genomeDir 041_d_STAR_mouseMT_reference \
     --genomeFastaFiles $G_FASTA \
     --sjdbGTFfile $G_GTF \
     --runThreadN 4 \
     --genomeSAindexNbases 5 \
     --sjdbOverhang 99 