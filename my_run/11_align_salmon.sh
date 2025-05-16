#!/usr/bin/bash
# pseudo alignment of the mouseMT reads with salmon


dataDIR=fastq_files

sourceFILE=sampleNames.txt


genomeDIR=MusMT_index

for fastqFILE in `cat $sourceFILE`
do
outDIR=033_d_salmon_mouseMT_${fastqFILE%.*}

 mkdir -p $outDIR

 salmon quant -i $genomeDIR -l A \
              -r $dataDIR/${fastqFILE}.fastq \
              -p 4 --validateMappings --gcBias --seqBias \
              -o $outDIR
done
