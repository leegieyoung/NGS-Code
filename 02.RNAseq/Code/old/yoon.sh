#!/bin/bash

for i in `ls`

do

mkdir $Colon/tmp/$i

hisat2 -q -p 68 -x $Ref/Hindexes/fusion_Homo_sapiens.GRCh38.79 \
        -1 $Raw/$i/$i'.R1.fq.gz' \
        -2 $Raw/$i/$i'.R2.fq.gz' \
        -S $Colon/tmp/$i/hisat2.dta.sam \
        --downstream-transcriptome-assembly
#        2> $Colon/hisat2/$i/hisat2.dta.log

done

samtools view -@ 68 -b hisat2.dta.sam | samtools sort -@ 68 -o hisat2.dta.sorted.bam

samtools index hisat2.dta.sorted.bam

cd /scratch/x1997a07/data/Raw_KR_kribb/IBD/high_QC

featureCounts -Q 20 -T 64 -a /scratch/x1997a07/data/reference/Homo_sapiens.GRCh38.79.gtf  -g 'gene_name' -o /scratch/x1997a07/data/RNA-seq/IBD/featurecount/IBD_ALL_symbol.txt -M *.bam
