#!/bin/sh
if [ $# -ne 4 ];then
        echo "Please enter merge_result/input1/input2/input3"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
OUTPUT=$1
Input1=$2
Input2=$3
Input3=$4
GWAS_path="/scratch/hpc46a05/GWAS/result"
code_path="/scratch/hpc46a05/GWAS/Code/plink_code_folder"
mkdir ${GWAS_path}/QC_Imputed_${OUTPUT}
result_folder="${GWAS_path}/QC_Imputed_${OUTPUT}"
#imputed_sample folder
raw_folder="/data/keeyoung/GWAS/raw_data/"
#beagle_folder="/scratch/x1997a11/GWAS/pdxen_AD/beagle/Result_folder"

#1.SampleQC(geno mind impute-sex hwe)
mkdir ${result_folder}/QC
mkdir ${result_folder}/MaMi


#.mergelist
Input1_result_folder="${GWAS_path}/QC_Imputed_${Input1}"
Input2_result_folder="${GWAS_path}/QC_Imputed_${Input2}"
Input3_result_folder="${GWAS_path}/QC_Imputed_${Input3}"

echo ${Input1_result_folder}/MaMi/QC_Imputed_${Input1}.bed ${Input1_result_folder}/MaMi/QC_Imputed_${Input1}.bim ${Input1_result_folder}/MaMi/QC_Imputed_${Input1}.fam > ${Input1_result_folder}/MaMi/QC_Imputed_${Input1}.mergelist
echo ${Input2_result_folder}/MaMi/QC_Imputed_${Input2}.bed ${Input2_result_folder}/MaMi/QC_Imputed_${Input2}.bim ${Input2_result_folder}/MaMi/QC_Imputed_${Input2}.fam > ${Input2_result_folder}/MaMi/QC_Imputed_${Input2}.mergelist
echo ${Input3_result_folder}/MaMi/QC_Imputed_${Input3}.bed ${Input3_result_folder}/MaMi/QC_Imputed_${Input3}.bim ${Input3_result_folder}/MaMi/QC_Imputed_${Input3}.fam > ${Input3_result_folder}/MaMi/QC_Imputed_${Input3}.mergelist

cat ${Input2_result_folder}/MaMi/QC_Imputed_${Input2}.mergelist ${Input3_result_folder}/MaMi/QC_Imputed_${Input3}.mergelist > ${result_folder}/QC/mergelist.txt

##missnp 가 없다면 아래 코드는 정상작동함
plink --bfile ${Input1_result_folder}/MaMi/QC_Imputed_${Input1} \
 --merge-list ${result_folder}/QC/mergelist.txt \
 --make-bed \
 --keep-allele-order \
 --out ${result_folder}/MaMi/QC_Imputed_${OUTPUT}


##missnp 가 있는 경우에만 작동함
plink --bfile ${Input1_result_folder}/MaMi/QC_Imputed_${Input1} \
 --exclude ${result_folder}/MaMi/QC_Imputed_${OUTPUT}-merge.missnp \
 --make-bed \
 --keep-allele-order \
 --out ${result_folder}/MaMi/QC_rmmissnp_Imputed_${Input1}

plink --bfile ${Input2_result_folder}/MaMi/QC_Imputed_${Input2} \
 --exclude ${result_folder}/MaMi/QC_Imputed_${OUTPUT}-merge.missnp \
 --make-bed \
 --keep-allele-order \
 --out ${result_folder}/MaMi/QC_rmmissnp_Imputed_${Input2}

plink --bfile ${Input3_result_folder}/MaMi/QC_Imputed_${Input3} \
 --exclude ${result_folder}/MaMi/QC_Imputed_${OUTPUT}-merge.missnp \
 --make-bed \
 --keep-allele-order \
 --out ${result_folder}/MaMi/QC_rmmissnp_Imputed_${Input3}

echo ${result_folder}/MaMi/QC_rmmissnp_Imputed_${Input1}.bed ${result_folder}/MaMi/QC_rmmissnp_Imputed_${Input1}.bim ${result_folder}/MaMi/QC_rmmissnp_Imputed_${Input1}.fam > ${result_folder}/MaMi/QC_rmmissnp_Imputed_${Input1}.mergelist
echo ${result_folder}/MaMi/QC_rmmissnp_Imputed_${Input2}.bed ${result_folder}/MaMi/QC_rmmissnp_Imputed_${Input2}.bim ${result_folder}/MaMi/QC_rmmissnp_Imputed_${Input2}.fam > ${result_folder}/MaMi/QC_rmmissnp_Imputed_${Input2}.mergelist
echo ${result_folder}/MaMi/QC_rmmissnp_Imputed_${Input3}.bed ${result_folder}/MaMi/QC_rmmissnp_Imputed_${Input3}.bim ${result_folder}/MaMi/QC_rmmissnp_Imputed_${Input3}.fam > ${result_folder}/MaMi/QC_rmmissnp_Imputed_${Input3}.mergelist

cat ${result_folder}/MaMi/QC_rmmissnp_Imputed_${Input2}.mergelist ${result_folder}/MaMi/QC_rmmissnp_Imputed_${Input3}.mergelist > ${result_folder}/QC/mergelist_rmmissnp.txt

plink --bfile ${result_folder}/MaMi/QC_rmmissnp_Imputed_${Input1} \
 --merge-list ${result_folder}/QC/mergelist_rmmissnp.txt \
 --exclude ${result_folder}/MaMi/QC_Imputed_${OUTPUT}-merge.missnp \
 --make-bed \
 --keep-allele-order \
 --out ${result_folder}/MaMi/QC_Imputed_${OUTPUT}

#=====================================merge End================================================

#=======================================rm NA =================================================
plink --bfile ${result_folder}/MaMi/QC_Imputed_${OUTPUT} \
 --freq case-control \
 --keep-allele-order \
 --out ${result_folder}/MaMi/NoQC_rmmissnp_Imputed_${OUTPUT}_freq


#After merge, and QC start

awk '{print $2, $5, $6}' ${result_folder}/MaMi/NoQC_rmmissnp_Imputed_${OUTPUT}_freq.frq.cc > ${result_folder}/MaMi/NonCheck_NA.snp
grep -v -w "NA" ${result_folder}/MaMi/NonCheck_NA.snp | awk '{print $1}' | sed '1d'  > ${result_folder}/MaMi/Check_NA.snp

plink --bfile ${result_folder}/MaMi/QC_Imputed_${OUTPUT} \
 --extract ${result_folder}/MaMi/Check_NA.snp \
 --keep-allele-order \
 --make-bed \
 --out ${result_folder}/MaMi/NoQC_rmIndep_rmmissnp_Imputed_${OUTPUT}

plink --bfile ${result_folder}/MaMi/NoQC_rmIndep_rmmissnp_Imputed_${OUTPUT} \
 --geno 0.05 \
 --keep-allele-order \
 --make-bed \
 --out ${result_folder}/MaMi/NoQC_g

plink --bfile ${result_folder}/MaMi/NoQC_g \
 --mind 0.2 \
 --keep-allele-order \
 --make-bed \
 --out ${result_folder}/MaMi/NoQC_g_m

plink --bfile ${result_folder}/MaMi/NoQC_g_m \
 --maf 0.05 \
 --keep-allele-order \
 --make-bed \
 --out ${result_folder}/MaMi/NoQC_g_m_maf

plink --bfile ${result_folder}/MaMi/NoQC_g_m_maf \
 --hwe 1e-6 \
 --keep-allele-order \
 --make-bed \
 --out ${result_folder}/MaMi/QC_Imputed_${OUTPUT}

echo ""
echo "=====================Delete Used Data==============="
echo ""
rm ${result_folder}/MaMi/NoQC*.bed
rm ${result_folder}/MaMi/NoQC*.bim
rm ${result_folder}/MaMi/NoQC*.fam
rm ${result_folder}/MaMi/QC_rmmissnp_Imputed_*.bed
rm ${result_folder}/MaMi/QC_rmmissnp_Imputed_*.bim
rm ${result_folder}/MaMi/QC_rmmissnp_Imputed_*.fam
rm ${result_folder}/MaMi/QC_rmmissnp_Imputed_*.mergelist
rm ${result_folder}/MaMi/*.nosex
rm ${result_folder}/MaMi/*.frq.cc
rm ${result_folder}/MaMi/*.snp
