##!/bin/sh
if [ $# -ne 3 ];then
        echo "Please enter Sample_Name"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
#beagle 진행시Same position 제거를 안했을 경우 사용하는 코드
Sample=$1
Thread=$2
MEM=$3
Dir="/ichrogene/project/temp/gylee/0.GWAS"
QC_dir="${Dir}/1.Input/${Sample}/1.QC"
MaMi_dir="${Dir}/1.Input/${Sample}/2.MaMi"
p2_Dir="/ichrogene/project/temp/gylee/Singurality/plink2"
Inversion="/ichrogene/project/temp/gylee/0.GWAS/REFERENCE/inversion.txt"
Output_dir="${Dir}/2.plink_result/${Sample}"

mkdir -p ${QC_dir}
mkdir -p ${Output_dir} 

export QC_DIR="${QC_dir}/"
export RESULT_DIR="${Output_dir}/"
echo ${ANA_DIR}
export SAMPLE="${Sample}"


for A in $(seq 1 22)
do
${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe \
 --memory ${MEM} \
 --threads ${Thread} \
 --make-bed \
 --chr ${A} \
 --out ${Output_dir}/chr${A}.imputation_g_m_maf_hwe_bfile 
done
