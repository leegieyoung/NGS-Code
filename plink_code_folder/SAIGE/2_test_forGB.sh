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
Output_dir="${Dir}/2.plink_result/${Sample}"

mkdir -p ${QC_dir}
mkdir -p ${Output_dir} 

export QC_DIR="${QC_dir}/"
export RESULT_DIR="${Output_dir}/"
#echo ${ANA_DIR}
export SAMPLE="${Sample}"


${p2_Dir}/plink2 --bfile ${QC_dir}/${Sample}_g_m_maf_hwe_het_KC1st \
 --glm sex hide-covar \
 --covar ${QC_dir}/covariate.txt \
 --covar-name PC1-PC10 \
 --covar-variance-standardize \
 --out ${Output_dir}/${Sample}

awk '{if($18 < 8.6e-8) print $0}' ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/raw_${Sample}.PHENO1.glm.logistic.hybrid.FDR
singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/Manhattan_plot.R

awk '{if($18 < 1e-5) print $0}' ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/raw_${Sample}.PHENO1.glm.logistic.hybrid.suggest

head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
cat ${Output_dir}/header ${Output_dir}/raw_${Sample}.PHENO1.glm.logistic.hybrid.FDR > ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid.FDR
cat ${Output_dir}/header ${Output_dir}/raw_${Sample}.PHENO1.glm.logistic.hybrid.suggest > ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid.suggest

rm -rf ${Output_dir}/raw_*

