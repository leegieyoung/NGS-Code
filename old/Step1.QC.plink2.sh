#!/bin/bash
Dir="/ichrogene/project/temp/gylee/0.GWAS/05.Imputation/0.INPUT"
p2_Dir="/ichrogene/project/temp/gylee/Singularity/plink2"

# 1. Remove duplicate
${p2_Dir}/plink2 --bfile /ichrogene/project/UKB_Chip/Ukb \
 --make-just-bim \
 --rm-dup \
 --out ${Dir}/rmdupli

# 2. Remove INDEL
awk '$5 != "-" && $6 != "-" && $5 != "I" && $6 != "I" && $5 != "D" && $6 != "D" && $5 != "0" && $6 != "0" && $5 != "N" && $6 != "N" && $5 != "." && $6 != "." ' ${Dir}/rmdupli.bim > ${Dir}/rmdupli_rmINDEL.bim

awk '{print $2}' ${Dir}/rmdupli_rmINDEL.bim > ${Dir}/rmdupli_rmINDEL.list

# 3. Only acgt
${p2_Dir}/plink2 --bfile /ichrogene/project/UKB_Chip/Ukb \
 --keep-allele-order \
 --make-pgen \
 --snps-only just-acgt \
 --extract ${Dir}/rmdupli_rmINDEL.list \
 --out ${Dir}/Ukb_rmdupli_rmINDEL

# 4. MAF and Keep, Keep 에 사용되는 샘플 리스트에 따라 아래 옵션 변경 필요
for Chr in $(seq 1 22)
do
${p2_Dir}/plink2 --pfile ${Dir}/Ukb/Ukb_rmdupli_rmINDEL \
 --recode vcf-iid bgz \
 --maf 0.01 \
 --chr ${Chr} \
 --keep /ichrogene/project/temp/gylee/0.GWAS/05.Imputation/0.INPUT/Cancers_67630/impulation_list.list \
 --out ${Dir}/Cancers_67630/Cancers_67630_${Chr}_rmdupli_rmINDEL
done
