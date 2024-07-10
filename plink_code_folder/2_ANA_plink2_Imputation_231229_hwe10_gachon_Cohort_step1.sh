#!/bin/sh
if [ $# -ne 4 ];then
        echo "Please enter Sample_Name"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
#beagle 진행시Same position 제거를 안했을 경우 사용하는 코드
Sample=$1
Thread=$2
MEM=$3
Cohort=$4
Dir="/ichrogene/project/temp/gylee/0.GWAS"
p2_Dir="/ichrogene/project/temp/gylee/Singularity/plink2"
Inversion="/ichrogene/project/temp/gylee/0.GWAS/REFERENCE/inversion.txt"
Output_dir="${Dir}/2.plink_result/${Sample}_${Cohort}"

mkdir -p ${Output_dir} 

export RESULT_DIR="${Output_dir}/"
#echo ${ANA_DIR}
export SAMPLE="${Sample}"

#cp /ichrogene/project/temp/gylee/0.GWAS/05.Imputation/2.minimac4/minimac4.1.0/${Sample}/${Sample}_pheno.txt ${Output_dir}/${Sample}_pheno.txt

cp /ichrogene/project/temp/gylee/0.GWAS/05.Imputation/2.minimac4/minimac4.1.0/${Sample}/${Sample}_sex.txt ${Output_dir}/${Sample}_sex.txt


${p2_Dir}/plink2 --pfile /ichrogene/project/temp/gylee/0.GWAS/05.Imputation/2.minimac4/minimac4.1.0/${Sample}/imputation \
 --keep ${Output_dir}/${Sample}_${Cohort}.list \
 --make-pgen \
 --rm-dup force-first \
 --memory ${MEM} \
 --update-sex ${Output_dir}/${Sample}_sex.txt \
 --out ${Output_dir}/imputation

${p2_Dir}/plink2 --pfile ${Output_dir}/imputation \
 --geno 0.02 \
 --make-pgen \
 --memory ${MEM} \
 --out ${Output_dir}/imputation_g

${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g \
 --mind 0.02 \
 --make-pgen \
 --memory ${MEM} \
 --out ${Output_dir}/imputation_g_m

${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m \
 --maf 0.01 \
 --make-pgen \
 --memory ${MEM} \
 --out ${Output_dir}/imputation_g_m_maf

${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf \
 --hwe 1e-50 \
 --memory ${MEM} \
 --make-pgen \
 --out ${Output_dir}/imputation_g_m_maf_hwe

${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe \
 --exclude ${Inversion} \
 --indep-pairwise 50 5 0.02 \
 --pheno ${Output_dir}/${Sample}_pheno.txt \
 --threads ${Thread} \
 --memory ${MEM} \
 --out ${Output_dir}/${Sample}_indepSNP

#het-missing.R
${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe \
 --missing \
 --out ${Output_dir}/Before_miss

${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe \
 --extract ${Output_dir}/${Sample}_indepSNP.prune.in \
 --het \
 --out ${Output_dir}/Before_het

awk '{print $4"\tPheno"}' ${Output_dir}/Before_miss.smiss > ${Output_dir}/B_F_MISS

paste -d '\t' ${Output_dir}/Before_het.het ${Output_dir}/B_F_MISS > ${Output_dir}/R_BeforeQC_check.txt

singularity exec --bind /ichrogene/:/ichrogene/ /ichrogene/project/temp/gylee/Singularity/GWAS.sif Rscript /ichrogene/project/temp/gylee/0.GWAS/Code/R/het-missing.R

#analysis
${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe \
 --double-id \
 --extract ${Output_dir}/${Sample}_indepSNP.prune.in \
 --pca \
 --memory ${MEM} \
 --set-missing-var-ids @:# \
 --threads ${Thread} \
 --out ${Output_dir}/${Sample}_PCA

singularity exec --bind /ichrogene/:/ichrogene/ /ichrogene/project/temp/gylee/Singularity/GWAS.sif Rscript /ichrogene/project/temp/gylee/0.GWAS/Code/R/covariate_PCA.R

#hand-work
singularity exec --bind /ichrogene/:/ichrogene/ /ichrogene/project/temp/gylee/Singularity/PRS.sif Rscript /ichrogene/project/temp/gylee/0.GWAS/Code/R/PCA_plot.R

##rm outlier
#awk '{print $1}' ${Output_dir}/covariate_forplink_PCA_outlier.txt > ${Output_dir}/covariate_forplink_PCA_outlier.list
#sed -i '1,1d' ${Output_dir}/covariate_forplink_PCA_outlier.list
#${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe \
# --keep ${Output_dir}/covariate_forplink_PCA_outlier.list \
# --make-pgen \
# --out ${Output_dir}/imputation_g_m_maf_hwe_rmoutlier
#
##het-missing.R
#${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe_rmoutlier \
# --missing \
# --out ${Output_dir}/Before_rmoutlier_miss
#
#${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe_rmoutlier \
# --extract ${Output_dir}/${Sample}_indepSNP.prune.in \
# --het \
# --out ${Output_dir}/Before_rmoutlier_het
#
#awk '{print $4"\tPheno"}' ${Output_dir}/Before_rmoutlier_miss.smiss > ${Output_dir}/B_F_MISS
#
#paste -d '\t' ${Output_dir}/Before_rmoutlier_het.het ${Output_dir}/B_F_MISS > ${Output_dir}/R_AfterQC_check.txt
#
#singularity exec --bind /ichrogene/:/ichrogene/ /ichrogene/project/temp/gylee/Singularity/GWAS.sif Rscript /ichrogene/project/temp/gylee/0.GWAS/Code/R/het-missing_AfterQC.R
#
#${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe_rmoutlier \
# --glm sex hide-covar \
# --memory ${MEM} \
# --pheno ${Output_dir}/${Sample}_pheno.txt \
# --covar ${Output_dir}/covariate_forplink_PCA_outlier.txt \
# --covar-name PC1-PC10 \
# --covar-variance-standardize \
# --out ${Output_dir}/${Sample}
#
#cp ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid ${Output_dir}/${Sample}.glm.repli
#cp ${Output_dir}/${Sample}.Pheno.glm.logistic.hybrid ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid
#
#grep -v "^#" ${Output_dir}/imputation_g_m_maf_hwe.pvar | awk '{print $3}' > ${Output_dir}/variant.info
#
##Bonferroni correction & Make FDR file
#singularity exec --bind /ichrogene/:/ichrogene/ /ichrogene/project/temp/gylee/Singularity/PRS.sif Rscript /ichrogene/project/temp/gylee/0.GWAS/Code/R/Manhattan_plot.R
