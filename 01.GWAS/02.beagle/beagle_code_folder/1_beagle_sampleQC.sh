#!/bin/bash
s1=$1
s2=$2
INPUT="/scratch/hpc46a05/GWAS/raw_data/${s1}"
OUTPUT="/scratch/hpc46a05/GWAS/beagle_result"

mkdir ${OUTPUT}/1.${s1}
#QC maf
plink --bfile ${INPUT}/${s1} --maf 0.001 --allow-no-sex --keep-allele-order --make-bed --out ${OUTPUT}/1.${s1}/${s2}_maf
plink --bfile ${OUTPUT}/1.${s1}/${s2}_maf --allow-no-sex --keep-allele-order --recode vcf-iid bgz --out ${OUTPUT}/1.${s1}/${s2}_maf

#Remove indels
zcat ${OUTPUT}/1.${s1}/${s2}_maf.vcf.gz | awk '$4 != "-" && $5 != "-" && $4 != "I" && $5 != "I" && $4 != "D" && $5 != "D" && $4 != "0" && $5 != "0" && $4 != "N" && $5 != "N" && $4 != "." && $5 != "." '  | gzip > ${OUTPUT}/1.${s1}/${s2}_maf_noINDEL.vcf.gz

#Remove duplicated SNPs
bcftools norm -d all -O z --threads 64 -o ${OUTPUT}/1.${s1}/${s2}_maf_noINDEL_rmDupID.vcf.gz ${OUTPUT}/1.${s1}/${s2}_maf_noINDEL.vcf.gz 

plink --vcf ${OUTPUT}/1.${s1}/${s2}_maf_noINDEL_rmDupID.vcf.gz \
 --double-id \
 --keep-allele-order \
 --allow-no-sex \
 --make-bed \
 --out ${OUTPUT}/1.${s1}/${s2}_maf_noINDEL_rmDupID_noliftover

#Extract Autosomal
grep -v '^26' ${OUTPUT}/1.${s1}/${s2}_maf_noINDEL_rmDupID_noliftover.bim > ${OUTPUT}/1.${s1}/rmother1.bim
grep -v '^25' ${OUTPUT}/1.${s1}/rmother1.bim > ${OUTPUT}/1.${s1}/rmother2.bim
grep -v '^24' ${OUTPUT}/1.${s1}/rmother2.bim > ${OUTPUT}/1.${s1}/rmother3.bim
grep -v '^23' ${OUTPUT}/1.${s1}/rmother3.bim > ${OUTPUT}/1.${s1}/rmother4.bim
grep -v '^MT' ${OUTPUT}/1.${s1}/rmother4.bim > ${OUTPUT}/1.${s1}/rmother5.bim
grep -v '^Y' ${OUTPUT}/1.${s1}/rmother5.bim > ${OUTPUT}/1.${s1}/rmother6.bim
grep -v '^X' ${OUTPUT}/1.${s1}/rmother6.bim > ${OUTPUT}/1.${s1}/${s2}_maf_noINDEL_rmDupID_noliftover.autosomal.bim

awk '{print $2}' ${OUTPUT}/1.${s1}/${s2}_maf_noINDEL_rmDupID_noliftover.autosomal.bim > ${OUTPUT}/1.${s1}/${s2}_maf_noINDEL_rmDupID_noliftover.autosomal.txt

plink --bfile ${OUTPUT}/1.${s1}/${s2}_maf_noINDEL_rmDupID_noliftover \
 --extract ${OUTPUT}/1.${s1}/${s2}_maf_noINDEL_rmDupID_noliftover.autosomal.txt \
 --allow-no-sex \
 --make-bed \
 --out ${OUTPUT}/1.${s1}/${s2}_maf_noINDEL_rmDupID_noliftover.autosomal

#Split Chromosomes
for chr in $(seq 1 22)
do
plink --bfile ${OUTPUT}/1.${s1}/${s2}_maf_noINDEL_rmDupID_noliftover.autosomal --allow-no-sex --keep-allele-order --chr $chr --double-id --recode vcf-iid bgz --out ${OUTPUT}/1.${s1}/CHR${chr}_${s2}_noIndel_rmDupID
bcftools index --threads 64 ${OUTPUT}/1.${s1}/CHR${chr}_${s2}_noIndel_rmDupID.vcf.gz
done
