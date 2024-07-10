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

#1.SampleQC(geno mind impute-sex hwe)
mkdir -p ${QC_dir}
mkdir -p ${MaMi_dir}

plink --bfile ${QC_dir}/${Sample} \
 --geno 0.02 \
 --make-bed \
 --threads ${Thread} \
 --out ${QC_dir}/${Sample}_g


plink --bfile ${QC_dir}/${Sample}_g \
 --mind 0.02 \
 --make-bed \
 --threads ${Thread} \
 --out ${QC_dir}/${Sample}_g_m

plink --bfile ${QC_dir}/${Sample}_g_m \
 --maf 0.01 \
 --make-bed \
 --threads ${Thread} \
 --out ${QC_dir}/${Sample}_g_m_maf

plink --bfile ${QC_dir}/${Sample}_g_m_maf \
 --hwe 1e-6 \
 --make-bed \
 --threads ${Thread} \
 --out ${QC_dir}/${Sample}_g_m_maf_hwe

echo ""
echo "===================Delete Used Data========================"
echo ""
