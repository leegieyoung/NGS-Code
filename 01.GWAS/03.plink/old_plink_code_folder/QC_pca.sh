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
Case_pheno="${analysis_folder}/merge/Case_pheno.txt"
Control_pheno="${analysis_folder}/merge/Control_pheno.txt"
inversion="/scratch/hpc46a05/REFERENCE/inversion.txt"
Rcode="/scratch/hpc46a05/GWAS/Rcode"

#PCA
mkdir -p ${QC_folder}/PCA
plink --bfile ${GWAS_path}/QC_${Sample}/MaMi/QC_${Sample} \
 --double-id \
 --pca 10 \
 --set-missing-var-ids @:# \
 --keep-allele-order \
 --out ${QC_folder}/PCA/QC_${Sample}_PCA

awk '{print $1, $3, $4, $5, $6}' ${QC_folder}/PCA/QC_${Sample}_PCA.eigenvec > ${QC_folder}/PCA/raw_QC_${Sample}_PCA.csv
sed -i '1i\name PC1 PC2 PC3 PC4' ${QC_folder}/PCA/raw_QC_${Sample}_PCA.csv
paste -d '\t' ${QC_folder}/PCA/raw_QC_${Sample}_PCA.csv ${QC_folder}/PCA/QC_${Sample}_PCA.eigenval > ${QC_folder}/PCA/raw_QC_${Sample}_PCA.eigenval.csv

awk '{print $6}' ${GWAS_path}/QC_${Sample}/MaMi/QC_${Sample}.fam > ${QC_folder}/PCA/QC_${Sample}.pheno
sed -i -e 's/2/Case/g' -e 's/1/Control/g' ${QC_folder}/PCA/QC_${Sample}.pheno
sed -i '1iname' ${QC_folder}/PCA/QC_${Sample}.pheno
awk '{print $1}' ${GWAS_path}/QC_${Sample}/MaMi/QC_${Sample}.fam > ${QC_folder}/PCA/QC_${Sample}.sample
sed -i '1isample' ${QC_folder}/PCA/QC_${Sample}.sample
awk '{print $3, $4}' ${QC_folder}/PCA/QC_${Sample}_PCA.eigenvec > ${QC_folder}/PCA/raw_QC_${Sample}_PCA.PC12
sed -i '1iPC1 PC2' ${QC_folder}/PCA/raw_QC_${Sample}_PCA.PC12
paste -d ' ' ${QC_folder}/PCA/QC_${Sample}.pheno ${QC_folder}/PCA/raw_QC_${Sample}_PCA.PC12 ${QC_folder}/PCA/QC_${Sample}.sample > ${QC_folder}/PCA/PCA.txt

cp ${Rcode}/pca.R ${QC_folder}/PCA/

