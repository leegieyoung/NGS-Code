#!/bin/sh
if [ $# -ne 3 ];then
        echo "Please enter Sample_Name"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
#beagle 진행시Same position 제거를 안했을 경우 사용하는 코드
Sample=$1
Thread=$2
MEM=$3
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

#cp /mnt/nas/gylee/0.GWAS/05.Imputation/2.minimac4/minimac4.1.0/${Sample}/${Sample}_pheno.txt ${Output_dir}/${Sample}_pheno.txt
#
#${p2_Dir}/plink2 --pfile /mnt/nas/gylee/0.GWAS/05.Imputation/2.minimac4/minimac4.1.0/${Sample}/imputation \
# --geno 0.02 \
# --make-pgen \
# --memory ${MEM} \
# --update-sex /mnt/nas/gylee/0.GWAS/1.Input/Ukb_obesity/1.raw/subset_of_ukb45411_for_sex.txt \
# --out ${Output_dir}/imputation_g
#
#${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g \
# --mind 0.02 \
# --make-pgen \
# --memory ${MEM} \
# --out ${Output_dir}/imputation_g_m
#
#${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m \
# --maf 0.01 \
# --make-pgen \
# --memory ${MEM} \
# --out ${Output_dir}/imputation_g_m_maf
#
#${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf \
# --hwe 1e-50 \
# --memory ${MEM} \
# --make-pgen \
# --out ${Output_dir}/imputation_g_m_maf_hwe
#
#${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe \
# --exclude ${Inversion} \
# --indep-pairwise 50 5 0.02 \
# --pheno ${Output_dir}/${Sample}_pheno.txt \
# --threads ${Thread} \
# --memory ${MEM} \
# --out ${Output_dir}/${Sample}_indepSNP
#
#${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe \
# --double-id \
# --extract ${Output_dir}/${Sample}_indepSNP.prune.in \
# --pca approx 10 \
# --memory ${MEM} \
# --set-missing-var-ids @:# \
# --threads ${Thread} \
# --out ${Output_dir}/${Sample}_PCA
#
#singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/covariate.R
#
#${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe \
# --glm sex hide-covar \
# --memory ${MEM} \
# --pheno ${Output_dir}/${Sample}_pheno.txt \
# --covar ${Output_dir}/covariate.txt \
# --covar-name PC1-PC10 AGE \
# --covar-variance-standardize \
# --out ${Output_dir}/${Sample}
#
#cp ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid ${Output_dir}/${Sample}.glm.repli
#cp ${Output_dir}/${Sample}.Pheno.glm.logistic.hybrid ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid
#
#grep -v "^#" ${Output_dir}/imputation_g_m_maf_hwe.pvar | awk '{print $3}' > ${Output_dir}/variant.info
#
##Bonferroni correction & Make FDR file
#singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/Manhattan_plot.R
#
singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/3_prune.sh ${Sample} ${Thread} ${MEM}
#singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker.R

