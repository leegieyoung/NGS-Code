#!/bin/sh
if [ $# -ne 1 ];then
        echo "Please enter Sample_Name"
               exit
fi
Sample=$1
GWAS_path="/scratch/hpc46a05/GWAS/result"
code_path="/scratch/hpc46a05/GWAS/Code/plink_code_folder"
QC_folder="${GWAS_path}/QC_${Sample}"
#imputed_sample folder
beagle_folder="/scratch/hpc46a05/GWAS/beagle_result"
Sample_folder="${beagle_folder}/3.${Sample}_imputed_${Sample}"
#분석결과를 담는 파일
mkdir ${GWAS_path}/${Sample}_analysis_folder
analysis_folder="/${GWAS_path}/${Sample}_analysis_folder"
Case_pheno="${QC_folder}/Case_pheno.txt"
Control_pheno="${QC_folder}/Control_pheno.txt"
inversion="/scratch/hpc46a05/REFERENCE/inversion.txt"
Rcode="/scratch/hpc46a05/GWAS/Rcode"
#=================================================================

mkdir ${QC_folder}
#=========What is Code ? =======================================
#PCA
mkdir ${QC_folder}/PCA
plink --bfile ${QC_folder}/MaMi/QC_${Sample} \
 --double-id \
 --pca 10 \
 --set-missing-var-ids @:# \
 --keep-allele-order \
 --out ${QC_folder}/PCA/${Sample}_NoNA_PCA

awk '{print $1, $3, $4, $5, $6}' ${QC_folder}/PCA/${Sample}_NoNA_PCA.eigenvec > ${QC_folder}/PCA/raw_${Sample}_PCA.csv
sed -i '1i\name PC1 PC2 PC3 PC4' ${QC_folder}/PCA/raw_${Sample}_PCA.csv
paste -d '\t' ${QC_folder}/PCA/raw_${Sample}_PCA.csv ${QC_folder}/PCA/${Sample}_NoNA_PCA.eigenval > ${QC_folder}/PCA/raw_${Sample}_PCA-eigenval.csv

