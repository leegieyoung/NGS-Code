#!/bin/sh
if [ $# -ne 6 ];then
	echo "#usage: Please enter Case, Control, CHR, Start position and End position, Gene_name"
		exit
fi
Case=$1
Control=$2
s2=$3
Start=$4
number=15000
s3=$((Start-number))
End=$5
s4=$((End+number))
s5=$6

GWAS_path="/scratch/x1997a11/GWAS/pdxen_AD/result_folder"
code_path="/scratch/x1997a11/GWAS/pdxen_AD/Code_folder"
Case_folder="${GWAS_path}/QC_${Case}"
Control_folder="${GWAS_path}/QC_${Control}"
#분석결과를 담을 파일
mkdir ${GWAS_path}/${Case}_${Control}_analysis_folder
analysis_folder="/${GWAS_path}/${Case}_${Control}_analysis_folder"
Case_pheno="/scratch/x1997a11/GWAS/pdxen_AD/reference_folder/${Case}.txt"
Control_pheno="/scratch/x1997a11/GWAS/pdxen_AD/reference_folder/${Control}.txt"
inversion="/scratch/x1997a11/GWAS/pdxen_AD/reference_folder/inversion.txt"

mkdir ${analysis_folder}/merge/haploview/${s5}
grep -w "^${s2}" ${analysis_folder}/merge/${Case}_${Control}_NoNA.bim | awk '{if ($4 >= '${s3}'  && $4 <= '${s4}' ) print $2 }' > ${analysis_folder}/merge/haploview/${s5}/${s5}_SNP.txt


plink --bfile ${analysis_folder}/merge/${Case}_${Control}_NoNA \
 --extract ${analysis_folder}/merge/haploview/${s5}/${s5}_SNP.txt \
 --make-bed \
 --out ${analysis_folder}/merge/haploview/${s5}/${s5}


#haploview  용
awk '$5 == "A" || $5 == "T" || $5 == "C" || $5 == "G" {print $0}' ${analysis_folder}/merge/haploview/${s5}/${s5}.bim | awk '$6 == "A" || $6 == "T" || $6 == "C" || $6 == "G" {print $2}' > ${analysis_folder}/merge/haploview/${s5}/${s5}_ATCG.txt

plink --bfile ${analysis_folder}/merge/haploview/${s5}/${s5} \
 --extract ${analysis_folder}/merge/haploview/${s5}/${s5}_ATCG.txt \
 --recodeHV \
 --freq case-control  \
 --out ${analysis_folder}/merge/haploview/${s5}/${s5}_ATCG

plink --bfile ${analysis_folder}/merge/haploview/${s5}/${s5} \
 --extract ${analysis_folder}/merge/haploview/${s5}/${s5}_ATCG.txt \
 --logistic \
 --covar ${analysis_folder}/merge/logistic.convar_mds.txt \
 --sex \
 --keep-allele-order \
 --hide-covar \
 --out ${analysis_folder}/merge/haploview/${s5}/${s5}_ATCG

#maf 5% 미만 제거
awk '$5 >= 0.05 || $6 >= 0.05 {print $2}' ${analysis_folder}/merge/haploview/${s5}/${s5}_ATCG_freq.frq.cc > ${analysis_folder}/merge/haploview/${s5}/${s5}_ATCG_maf005.txt
plink --bfile ${analysis_folder}/merge/haploview/${s5}/${s5}_ATCG \
 --extract ${analysis_folder}/merge/haploview/${s5}/${s5}_ATCG_maf005.txt \
 --make-bed \
 --keep-allele-order \
 --out ${analysis_folder}/merge/haploview/${s5}/${s5}_ATCG_maf005

plink --bfile ${analysis_folder}/merge/haploview/${s5}/${s5}_ATCG_maf005 \
 --recodeHV \
 --keep-allele-order \
 --out ${analysis_folder}/merge/haploview/${s5}/${s5}_ATCG_maf005_recodeHV

plink --bfile ${analysis_folder}/merge/haploview/${s5}/${s5}_ATCG_maf005 \
 --freq case-control  \
 --keep-allele-order \
 --out ${analysis_folder}/merge/haploview/${s5}/${s5}_ATCG_maf005_freq

plink --bfile ${analysis_folder}/merge/haploview/${s5}/${s5}_ATCG_maf005 \
 --logistic \
 --covar ${analysis_folder}/merge/logistic.convar_mds.txt \
 --sex \
 --keep-allele-order \
 --hide-covar \
 --out ${analysis_folder}/merge/haploview/${s5}/${s5}_ATCG_maf005_logistic

awk '{print $2, $9}' ${analysis_folder}/merge/haploview/${s5}/${s5}_ATCG_maf005_logistic.assoc.logistic > ${analysis_folder}/merge/haploview/${s5}/${s5}_ATCG_maf005_rsID_logi
awk '{print $3, $4 , $5, $6, $7, $8}' ${analysis_folder}/merge/haploview/${s5}/${s5}_ATCG_maf005_freq.frq.cc > ${analysis_folder}/merge/haploview/${s5}/${s5}_ATCG_maf005_A1_A2_MAF_NCHROBS

paste ${analysis_folder}/merge/haploview/${s5}/${s5}_ATCG_maf005_rsID_logi ${analysis_folder}/merge/haploview/${s5}/${s5}_ATCG_maf005_A1_A2_MAF_NCHROBS > ${analysis_folder}/merge/haploview/${s5}/${s5}_ATCG_maf005_result.txt
sed -i 's/\t/ /g' ${analysis_folder}/merge/haploview/${s5}/${s5}_ATCG_maf005_result.txt
sed -i 's/  / /g' ${analysis_folder}/merge/haploview/${s5}/${s5}_ATCG_maf005_result.txt
