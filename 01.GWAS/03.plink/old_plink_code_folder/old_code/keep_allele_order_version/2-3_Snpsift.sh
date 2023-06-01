#!/bin/sh
if [ $# -ne 2 ];then
        echo "Please enter Sample_Name"
               exit
fi
Sample=$1
chr=$2
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
#snpsift dbsnp154 multi-core
java -Xmx64g -jar /scratch/x1997a11/SnpEff/snpEff/SnpSift.jar annotate /scratch/x1997a11/REFERENCE/Anno/dbsnp154/dbsnp154.vcf.gz ${analysis_folder}/merge/logistic/anno/raw_chr${chr}_rmsnp.vcf > ${analysis_folder}/merge/logistic/anno/chr${chr}_${Sample}_NoNA.vcf
