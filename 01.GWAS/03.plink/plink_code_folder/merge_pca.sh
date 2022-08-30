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
#=================================================================


awk '{print $1, $3, $4, $5, $6}' ${analysis_folder}/merge/PCA/${Sample}_NoNA_PCA.eigenvec > ${analysis_folder}/merge/PCA/raw_${Sample}_PCA.csv
sed -i '1i\name PC1 PC2 PC3 PC4' ${analysis_folder}/merge/PCA/raw_${Sample}_PCA.csv
paste -d '\t' ${analysis_folder}/merge/PCA/raw_${Sample}_PCA.csv ${analysis_folder}/merge/PCA/${Sample}_NoNA_PCA.eigenval > ${analysis_folder}/merge/PCA/raw_${Sample}_PCA-eigenval.csv

awk '{print $6}' ${analysis_folder}/merge/${Sample}_NoNA.fam > ${analysis_folder}/merge/PCA/${Sample}_NoNA.pheno
sed -i -e 's/2/Case/g' -e 's/1/Control/g' ${analysis_folder}/merge/PCA/${Sample}_NoNA.pheno
sed -i '1iname' ${analysis_folder}/merge/PCA/${Sample}_NoNA.pheno
awk '{print $1}' ${analysis_folder}/merge/${Sample}_NoNA.fam > ${analysis_folder}/merge/PCA/${Sample}_NoNA.sample
sed -i '1isample' ${analysis_folder}/merge/PCA/${Sample}_NoNA.sample
awk '{print $3, $4}' ${analysis_folder}/merge/PCA/${Sample}_NoNA_PCA.eigenvec > ${analysis_folder}/merge/PCA/raw_${Sample}_PCA.PC12
sed -i '1iPC1 PC2' ${analysis_folder}/merge/PCA/raw_${Sample}_PCA.PC12
paste -d ' ' ${analysis_folder}/merge/PCA/${Sample}_NoNA.pheno ${analysis_folder}/merge/PCA/raw_${Sample}_PCA.PC12 ${analysis_folder}/merge/PCA/${Sample}_NoNA.sample > ${analysis_folder}/merge/PCA/PCA.txt

cp ${Rcode}/pca.R ${analysis_folder}/merge/PCA/
