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

head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
cat ${Output_dir}/header ${Output_dir}/raw_${Sample}.PHENO1.glm.logistic.hybrid.FDR > ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid.FDR

awk '{print $3}' ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid.FDR > ${Output_dir}/${Sample}.snplist

singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/vep/vep.sif vep \
 -i ${Output_dir}/${Sample}.snplist --assembly GRCh37 --format id -o ${Output_dir}/${Sample}.snplist.vep --show_ref_allele --symbol --database --fork ${Thread}

singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/3_VEP_uniq_gene.R

mkdir ${Output_dir}/prune
${p2_Dir}/plink2 --bfile ${QC_dir}/${Sample}_g_m_maf_hwe_het_KC1st \
 --extract ${Output_dir}/${Sample}.snplist \
 --indep-pairwise 1000kb 0.2 \
 --out ${Output_dir}/prune/min0.2

rm -rf ${Output_dir}/prune/raw_prune.glm
touch ${Output_dir}/prune/raw_prune.glm

for A in $(cat ${Output_dir}/prune/min0.2.prune.in)
do
grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid.FDR >> ${Output_dir}/prune/raw_prune.glm
done
cat ${Output_dir}/header ${Output_dir}/prune/raw_prune.glm > ${Output_dir}/prune/prune.glm

singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/3_VEP_uniq_gene.R
