#!/bin/sh
if [ $# -ne 1 ];then
        echo "Please enter Input file name"
               exit
fi
Input=$1
Dir="/ichrogene/project/temp/gylee/0.GWAS/05.Imputation" #Change to your path
Image="/ichrogene/project/temp/gylee/Singularity"  #Change to your image file path
p2_Dir="/ichrogene/project/temp/gylee/Singularity/plink2"
scratch="/ichrogene"
ref_Dir="/ichrogene/project/temp/gylee/0.GWAS/REFERENCE/Imputation/reference"

mkdir ${Dir}/temp
mkdir ${Dir}/1.shapeit #Result folder of shapit4 that phasing tool
mkdir ${Dir}/2.minimac4  #Result folder of minimac4 that impuation tool

# 4. Divide by chromosome
echo "==============="
echo "Divide"
echo "=============="
for Chr in $(seq 1 22)
do
echo "==============="
echo "Chr${Chr} start"
echo "=============="

bcftools index -f --threads 32 ${Dir}/0.INPUT/${Input}/imputation_${Input}_${Chr}_rmdupli_rmINDEL.vcf.gz
done

