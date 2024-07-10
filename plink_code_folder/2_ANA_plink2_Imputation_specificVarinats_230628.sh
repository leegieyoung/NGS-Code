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

specific_dir="${Dir}/2.plink_result/Imputed_Ukb_T2D_ratio2"

#cp /mnt/nas/gylee/0.GWAS/05.Imputation/2.minimac4/minimac4.1.0/${Sample}/${Sample}_pheno.txt ${Output_dir}/${Sample}_pheno.txt
#cp ${specific_dir}/*_pheno.txt ${Output_dir}/${Sample}_pheno.txt
#awk '{print $1}' ${Output_dir}/${Sample}_pheno.txt > ${Output_dir}/${Sample}_sample.txt
#sed -i 's/IID/#IID/g' ${Output_dir}/${Sample}_sample.txt
#${p2_Dir}/plink2 --pfile /mnt/nas/gylee/0.GWAS/05.Imputation/2.minimac4/minimac4.1.0/Imputed_Ukb_obesity_BMI_WC/imputation \
# --make-pgen \
# --extract ${specific_dir}/specificVariants.list \
# --memory ${MEM} \
# --keep ${Output_dir}/${Sample}_sample.txt \
# --update-sex /mnt/nas/gylee/0.GWAS/1.Input/Ukb_obesity/1.raw/subset_of_ukb45411_for_sex.txt \
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

#Bonferroni correction & Make FDR file
#singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/Manhattan_plot.R

bash /mnt/nas/gylee/0.GWAS/Code/plink_code_folder/3_prune.sh Imputed_Ukb_obesity_BMI_WC_premium ${Thread} ${MEM}
bash /mnt/nas/gylee/0.GWAS/Code/plink_code_folder/4_patent_Obesity_T2D.sh ${Sample} ${Thread} ${MEM}
