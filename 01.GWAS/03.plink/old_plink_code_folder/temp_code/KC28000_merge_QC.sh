#!/bin/sh
if [ $# -ne 4 ];then
        echo "Please enter merge_result/input1/input2/input3"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
Sample=$1
Input1=$2
Input2=$3
Input3=$4
GWAS_path="/data/keeyoung/GWAS/result"
code_path="/scratch/x1997a11/GWAS/pdxen_AD/Code_folder"
mkdir ${GWAS_path}/QC_Imputed_${Sample}
result_folder="${GWAS_path}/QC_Imputed_${Sample}"
#imputed_sample folder
raw_folder="/data/keeyoung/GWAS/raw_data/"
#beagle_folder="/scratch/x1997a11/GWAS/pdxen_AD/beagle/Result_folder"

#1.SampleQC(geno mind impute-sex hwe)
mkdir ${result_folder}/raw_merge
mkdir ${result_folder}/QC
mkdir ${result_folder}/MaMi


#.mergelist
Input1_result_folder="${GWAS_path}/QC_Imputed_${Input1}"
Input2_result_folder="${GWAS_path}/QC_Imputed_${Input2}"
Input3_result_folder="${GWAS_path}/QC_Imputed_${Input3}"

echo ${Input1_result_folder}/QC/raw_NoQC_Imputed_${Input1}.bed ${Input1_result_folder}/QC/raw_NoQC_Imputed_${Input1}.bim ${Input1_result_folder}/QC/raw_NoQC_Imputed_${Input1}.fam > ${Input1_result_folder}/QC/raw_NoQC_Imputed_${Input1}.mergelist
echo ${Input2_result_folder}/QC/raw_NoQC_Imputed_${Input2}.bed ${Input2_result_folder}/QC/raw_NoQC_Imputed_${Input2}.bim ${Input2_result_folder}/QC/raw_NoQC_Imputed_${Input2}.fam > ${Input2_result_folder}/QC/raw_NoQC_Imputed_${Input2}.mergelist
echo ${Input3_result_folder}/QC/raw_NoQC_Imputed_${Input3}.bed ${Input3_result_folder}/QC/raw_NoQC_Imputed_${Input3}.bim ${Input3_result_folder}/QC/raw_NoQC_Imputed_${Input3}.fam > ${Input3_result_folder}/QC/raw_NoQC_Imputed_${Input3}.mergelist

cat ${Input2_result_folder}/QC/raw_NoQC_Imputed_${Input2}.mergelist ${Input3_result_folder}/QC/raw_NoQC_Imputed_${Input3}.mergelist > ${result_folder}/QC/mergelist.txt

plink --bfile ${Input1_result_folder}/QC/raw_NoQC_Imputed_${Input1} \
 --merge-list ${result_folder}/QC/mergelist.txt \
 --make-bed \
 --keep-allele-order \
 --out ${result_folder}/QC/raw_NoQC_Imputed_${Sample}

#missnp 가 없다면 아래 코드는 정상작동함
mv ${result_folder}/QC/raw_NoQC_Imputed_${Sample}.bed ${result_folder}/QC/raw_NoQC_rmmissnp_Imputed_${Sample}.bed
mv ${result_folder}/QC/raw_NoQC_Imputed_${Sample}.bim ${result_folder}/QC/raw_NoQC_rmmissnp_Imputed_${Sample}.bim
mv ${result_folder}/QC/raw_NoQC_Imputed_${Sample}.fam ${result_folder}/QC/raw_NoQC_rmmissnp_Imputed_${Sample}.fam

#missnp 가 있는 경우에만 작동함

plink --bfile ${Input1_result_folder}/QC/raw_NoQC_Imputed_${Input1} \
 --exclude ${result_folder}/QC/raw_NoQC_Imputed_${Sample}-merge.missnp \
 --make-bed \
 --keep-allele-order \
 --out ${result_folder}/QC/raw_NoQC_rmmissnp_Imputed_${Input1}

plink --bfile ${Input2_result_folder}/QC/raw_NoQC_Imputed_${Input2} \
 --exclude ${result_folder}/QC/raw_NoQC_Imputed_${Sample}-merge.missnp \
 --make-bed \
 --keep-allele-order \
 --out ${result_folder}/QC/raw_NoQC_rmmissnp_Imputed_${Input2}

plink --bfile ${Input3_result_folder}/QC/raw_NoQC_Imputed_${Input3} \
 --exclude ${result_folder}/QC/raw_NoQC_Imputed_${Sample}-merge.missnp \
 --make-bed \
 --keep-allele-order \
 --out ${result_folder}/QC/raw_NoQC_rmmissnp_Imputed_${Input3}

echo ${result_folder}/QC/raw_NoQC_rmmissnp_Imputed_${Input1}.bed ${result_folder}/QC/raw_NoQC_rmmissnp_Imputed_${Input1}.bim ${result_folder}/QC/raw_NoQC_rmmissnp_Imputed_${Input1}.fam > ${result_folder}/QC/raw_NoQC_rmmissnp_Imputed_${Input1}.mergelist
echo ${result_folder}/QC/raw_NoQC_rmmissnp_Imputed_${Input2}.bed ${result_folder}/QC/raw_NoQC_rmmissnp_Imputed_${Input2}.bim ${result_folder}/QC/raw_NoQC_rmmissnp_Imputed_${Input2}.fam > ${result_folder}/QC/raw_NoQC_rmmissnp_Imputed_${Input2}.mergelist
echo ${result_folder}/QC/raw_NoQC_rmmissnp_Imputed_${Input3}.bed ${result_folder}/QC/raw_NoQC_rmmissnp_Imputed_${Input3}.bim ${result_folder}/QC/raw_NoQC_rmmissnp_Imputed_${Input3}.fam > ${result_folder}/QC/raw_NoQC_rmmissnp_Imputed_${Input3}.mergelist

cat ${result_folder}/QC/raw_NoQC_rmmissnp_Imputed_${Input2}.mergelist ${result_folder}/QC/raw_NoQC_rmmissnp_Imputed_${Input3}.mergelist > ${result_folder}/QC/mergelist_rmmissnp.txt

plink --bfile ${result_folder}/QC/raw_NoQC_rmmissnp_Imputed_${Input1} \
 --merge-list ${result_folder}/QC/mergelist_rmmissnp.txt \
 --exclude ${result_folder}/QC/raw_NoQC_Imputed_${Sample}-merge.missnp \
 --make-bed \
 --keep-allele-order \
 --out ${result_folder}/QC/raw_NoQC_rmmissnp_Imputed_${Sample}

# QC
plink --bfile ${result_folder}/QC/raw_NoQC_rmmissnp_Imputed_${Sample} \
 --geno 0.2 \
 --keep-allele-order \
 --make-bed \
 --out ${result_folder}/QC/raw_NoQC_Imputed_${Sample}_g

plink --bfile ${result_folder}/QC/raw_NoQC_Imputed_${Sample}_g \
 --mind 0.2 \
 --keep-allele-order \
 --make-bed \
 --out ${result_folder}/QC/raw_NoQC_Imputed_${Sample}_g_m

plink --bfile ${result_folder}/QC/raw_NoQC_Imputed_${Sample}_g_m \
 --maf 0.01 \
 --keep-allele-order \
 --make-bed \
 --out ${result_folder}/QC/raw_NoQC_Imputed_${Sample}_g_m_maf

plink --bfile ${result_folder}/QC/raw_NoQC_Imputed_${Sample}_g_m_maf \
 --hwe 1e-6 \
 --keep-allele-order \
 --make-bed \
 --out ${result_folder}/QC/raw_NoQC_Imputed_${Sample}_g_m_maf_hwe

mv ${result_folder}/QC/raw_NoQC_Imputed_${Sample}_g_m_maf_hwe.bed ${result_folder}/MaMi/QC_Imputed_${Sample}.bed
mv ${result_folder}/QC/raw_NoQC_Imputed_${Sample}_g_m_maf_hwe.bim ${result_folder}/MaMi/QC_Imputed_${Sample}.bim
mv ${raw_folder}/${Sample}/${Sample}.fam ${result_folder}/MaMi/QC_Imputed_${Sample}.fam
