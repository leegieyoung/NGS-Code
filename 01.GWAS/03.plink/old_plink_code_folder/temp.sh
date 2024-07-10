#!/bin/sh
if [ $# -ne 3 ];then
        echo "Please enter merge_result/input1/input2/input3"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
OUTPUT=$1
Input1=$2
Input2=$3
GWAS_path="/scratch/hpc46a05/GWAS/result"
code_path="/scratch/hpc46a05/GWAS/Code/plink_code_folder"
mkdir ${GWAS_path}/QC_Imputed_${OUTPUT}
result_folder="${GWAS_path}/QC_Imputed_${OUTPUT}"
#imputed_sample folder
raw_folder="/data/keeyoung/GWAS/raw_data/"
#beagle_folder="/scratch/x1997a11/GWAS/pdxen_AD/beagle/Result_folder"



#=====================================merge End================================================

#=======================================rm NA =================================================
plink --bfile ${result_folder}/MaMi/raw_QC_Imputed_${OUTPUT} \
 --freq case-control \
 --keep-allele-order \
 --allow-no-sex \
 --out ${result_folder}/MaMi/NoQC_rmmissnp_Imputed_${OUTPUT}_freq


#After merge, and QC start

awk '{print $2, $5, $6}' ${result_folder}/MaMi/NoQC_rmmissnp_Imputed_${OUTPUT}_freq.frq.cc > ${result_folder}/MaMi/NonCheck_NA.snp
grep -v -w "NA" ${result_folder}/MaMi/NonCheck_NA.snp | awk '{print $1}' | sed '1d'  > ${result_folder}/MaMi/Check_NA.snp

plink --bfile ${result_folder}/MaMi/raw_QC_Imputed_${OUTPUT} \
 --extract ${result_folder}/MaMi/Check_NA.snp \
 --keep-allele-order \
 --make-bed \
 --allow-no-sex \
 --out ${result_folder}/MaMi/NoQC_rmIndep_rmmissnp_Imputed_${OUTPUT}

plink --bfile ${result_folder}/MaMi/NoQC_rmIndep_rmmissnp_Imputed_${OUTPUT} \
 --geno 0.2 \
 --allow-no-sex \
 --keep-allele-order \
 --make-bed \
 --out ${result_folder}/MaMi/NoQC_g

plink --bfile ${result_folder}/MaMi/NoQC_g \
 --mind 0.2 \
 --allow-no-sex \
 --keep-allele-order \
 --make-bed \
 --out ${result_folder}/MaMi/NoQC_g_m

plink --bfile ${result_folder}/MaMi/NoQC_g_m \
 --maf 0.05 \
 --allow-no-sex \
 --keep-allele-order \
 --make-bed \
 --out ${result_folder}/MaMi/NoQC_g_m_maf

plink --bfile ${result_folder}/MaMi/NoQC_g_m_maf \
 --hwe 1e-6 \
 --allow-no-sex \
 --keep-allele-order \
 --make-bed \
 --out ${result_folder}/MaMi/QC_Imputed_${OUTPUT}

echo ""
echo "=====================Delete Used Data==============="
echo ""
#rm ${result_folder}/MaMi/NoQC*.bed
#rm ${result_folder}/MaMi/NoQC*.bim
#rm ${result_folder}/MaMi/NoQC*.fam
#rm ${result_folder}/MaMi/QC_rmmissnp_Imputed_*.bed
#rm ${result_folder}/MaMi/QC_rmmissnp_Imputed_*.bim
#rm ${result_folder}/MaMi/QC_rmmissnp_Imputed_*.fam
#rm ${result_folder}/MaMi/QC_rmmissnp_Imputed_*.mergelist
#rm ${result_folder}/MaMi/*.nosex
#rm ${result_folder}/MaMi/*.frq.cc
#rm ${result_folder}/MaMi/*.snp
