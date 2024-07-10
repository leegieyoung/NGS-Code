#!/bin/sh
if [ $# -ne 1 ];then
        echo "Please enter Sample_Name"
               exit
fi
Sample=$1
cutoff=$2
GWAS_path="/scratch/x1997a11/GWAS/pdxen_AD/result_folder"
code_path="/scratch/x1997a11/GWAS/pdxen_AD/Code_folder"
Sample_folder="${GWAS_path}/QC_${Sample}"
#분석결과를 담을 파일
mkdir ${GWAS_path}/${Sample}_analysis_folder
analysis_folder="/${GWAS_path}/${Sample}_analysis_folder"
Case_pheno="${analysis_folder}/merge/Case_pheno.txt"
Control_pheno="${analysis_folder}/merge/Control_pheno.txt"
inversion="/scratch/x1997a11/GWAS/pdxen_AD/reference_folder/inversion.txt"
#=================================================================


#anno
mkdir ${analysis_folder}/merge/logistic/anno

plink --bfile ${analysis_folder}/merge/${Sample}_NoNA \
 --recode vcf-iid \
 --keep-allele-order \
 --double-id \
 --out ${analysis_folder}/merge/logistic/anno/raw_snp_vcf

awk '{print $1, $2, $3, $4, $5, $6, $7, $8}' ${analysis_folder}/merge/logistic/anno/raw_snp_vcf.vcf >  ${analysis_folder}/merge/logistic/anno/raw_snp_vcf_rmSam.vcf
grep '^#' ${analysis_folder}/merge/logistic/anno/raw_snp_vcf_rmSam.vcf > ${analysis_folder}/merge/logistic/anno/head.vcf
grep -v '^#' ${analysis_folder}/merge/logistic/anno/raw_snp_vcf_rmSam.vcf > ${analysis_folder}/merge/logistic/anno/nohead.vcf

awk '$3="." {print $0}' ${analysis_folder}/merge/logistic/anno/nohead.vcf > ${analysis_folder}/merge/logistic/anno/nohead_rmsnp.vcf
sed -i 's/ /\t/g' ${analysis_folder}/merge/logistic/anno/nohead_rmsnp.vcf
cat ${analysis_folder}/merge/logistic/anno/head.vcf ${analysis_folder}/merge/logistic/anno/nohead_rmsnp.vcf > ${analysis_folder}/merge/logistic/anno/rmsnp.vcf

for A in $(seq 1 22)
do
grep -w "^${A}" ${analysis_folder}/merge/logistic/anno/nohead_rmsnp.vcf > ${analysis_folder}/merge/logistic/anno/raw_chr${A}_nohead_rmsnp.vcf
cat ${analysis_folder}/merge/logistic/anno/head.vcf ${analysis_folder}/merge/logistic/anno/raw_chr${A}_nohead_rmsnp.vcf > ${analysis_folder}/merge/logistic/anno/raw_chr${A}_rmsnp.vcf
done

#snpsift dbsnp154

