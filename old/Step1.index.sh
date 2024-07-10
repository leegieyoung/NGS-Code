#!/bin/sh
Image="/ichrogene/project/temp/gylee/Singularity/"
Scratch="/ichrogene/"
Dir="/ichrogene/project/temp/gylee/0.GWAS/05.Imputation/0.INPUT/Ukb"
p2_Dir="/ichrogene/project/temp/gylee/Singularity/plink2"

# 4. Divide by chromosome
echo "==============="
echo "Divide"
echo "=============="
for Chr in $(seq 1 22)
do
echo "==============="
echo "Chr${Chr} start"
echo "=============="
singularity exec --bind ${Scratch}:${Scratch} ${Image}/minimac4.1.2.sif bcftools index --threads 32 -f ${Dir}/Ukb_${Chr}_rmdupli_rmINDEL.vcf.gz
done


