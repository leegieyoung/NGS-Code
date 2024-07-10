#!/bin/bash
if [ $# -ne 1 ];then
        echo "Please enter Input file name & Chromosme(Only number)"
               exit
fi
Input=$1
scratch="/ichrogene/"
Obj_Dir="/ichrogene/project/temp/gylee/REFERENCE"
Image="/ichrogene/project/temp/gylee/Singularity"  #Change to your image file path
Dir="/ichrogene/project/temp/gylee/0.GWAS/05.Imputation"
p2_Dir="/ichrogene/project/temp/gylee/Singularity/plink2"

for Chr in $(seq 1 22)
do
bcftools --version
nohup bcftools index -f ${Dir}/0.INPUT/${Input}/imputation_${Input}_${Chr}_rmdupli_rmINDEL.vcf.gz > /dev/null 2>&1 &
done

