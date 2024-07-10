#!/bin/sh
if [ $# -ne 2 ];then
        echo "Please enter Sample_Name"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
#beagle 진행시Same position 제거를 안했을 경우 사용하는 코드
Sample=$1
Thread=$2
Dir="/mnt/nas/gylee/0.GWAS"
QC_dir="${Dir}/1.Input/${Sample}/1.QC"
MaMi_dir="${Dir}/1.Input/${Sample}/2.MaMi"
p2_Dir="/mnt/nas/gylee/Singurality/plink2"
Inversion="/mnt/nas/gylee/0.GWAS/REFERENCE/inversion.txt"

#1.SampleQC(geno mind impute-sex hwe)
mkdir -p ${QC_dir}
mkdir -p ${MaMi_dir}

echo ${QC_dir}
##R pwd
export QC_DIR="${QC_dir}/"
echo ${ANA_DIR}
export SAMPLE="${Sample}"

#Remove 1st
${p2_Dir}/plink2 --bfile ${QC_dir}/${Sample}_g_m_maf_hwe_het \
 --king-cutoff 0.177 \
 --make-bed \
 --pheno ${QC_dir}/${Sample}_pheno.txt \
 --out ${QC_dir}/${Sample}_g_m_maf_hwe_het_KC1st
