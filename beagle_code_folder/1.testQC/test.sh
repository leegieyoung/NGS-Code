#!/bin/sh
INPUT=$1
#NAS="/mnt/nas"
NAS="/src/data"
dir="${NAS}/0.GWAS/Code/beagle_code_folder/1.testQC"

#QC_maf
#plink --file ${dir}/${INPUT} --maf 0.001 --make-bed --out ${dir}/${INPUT}_maf #ped, map
plink --bfile ${dir}/${INPUT} \
	--maf 0.001 --make-bed \
	--keep-allele-order \
	--allow-no-sex \
	--out ${dir}/${INPUT}_maf #bed bim fam

#Remove_indel
plink --bfile ${dir}/${INPUT}_maf \
	--recode vcf-iid bgz \
	--keep-allele-order \
	--allow-no-sex \
	--out ${dir}/${INPUT}_maf

zcat ${dir}/${INPUT}_maf.vcf.gz | awk '$4 != "-" && $5 != "-" && $4 != "I" && $5 != "I" && $4 != "D" && $5 != "D" && $4 != "0" && $5 != "0" && $4 != "N" && $5 != "N" && $4 != "." && $5 != "." ' | gzip > ${dir}/${INPUT}_maf_noINDEL.vcf.gz
bcftools norm -d all -O z \
	--threads 64 \
	-o ${dir}/${INPUT}_maf_noINDEL_rmDupID.vcf.gz ${dir}/${INPUT}_maf_noINDEL.vcf.gz

#Remove duplicated SNPs
plink --vcf ${dir}/${INPUT}_maf_noINDEL_rmDupID.vcf.gz \
	--double-id --make-bed \
	--keep-allele-order \
	--allow-no-sex \
	--out ${dir}/${INPUT}_maf_noINDEL_rmDupID

#Split Chromosomes
for chr in $(seq 1 22)
do
plink --bfile ${dir}/${INPUT}_maf_noINDEL_rmDupID \
	--chr ${chr} --recode vcf-iid bgz \
	--keep-allele-order \
	--allow-no-sex \
	--out ${dir}/CHR${chr}_${INPUT}_maf_noINDEL_rmDupID

bcftools index ${dir}/CHR${chr}_${INPUT}_maf_noINDEL_rmDupID.vcf.gz
done

