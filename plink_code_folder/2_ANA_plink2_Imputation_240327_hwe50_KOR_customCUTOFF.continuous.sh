##!/bin/sh
if [ $# -ne 4 ];then
        echo "Please enter Sample_Name"
               exit
fi
Sample=$1
Thread=$2
MEM=$3
TrainVali=$4
Dir="/ichrogene/project/temp/gylee/0.GWAS"
QC_dir="${Dir}/1.Input/${Sample}/1.QC"
MaMi_dir="${Dir}/1.Input/${Sample}/2.MaMi"
p2_Dir="/ichrogene/project/temp/gylee/Singularity/plink2"
Inversion="/ichrogene/project/temp/gylee/0.GWAS/REFERENCE/inversion.txt"
Output_dir="${Dir}/2.plink_result/${TrainVali}/${Sample}"

mkdir -p ${QC_dir}
mkdir -p ${Output_dir} 

export QC_DIR="${QC_dir}/"
export RESULT_DIR="${Output_dir}/"
echo ${ANA_DIR}
export SAMPLE="${Sample}"

#cp /ichrogene/project/temp/gylee/0.GWAS/05.Imputation/2.minimac4/minimac4.1.0/${Sample}/${Sample}_pheno.txt ${Output_dir}/${Sample}_pheno.txt

#${p2_Dir}/plink2 --pfile /ichrogene/project/temp/gylee/0.GWAS/05.Imputation/2.minimac4/minimac4.1.0/QC_Imputed_KNIH_72295ea_info0.8/QC_Imputed_KNIH_72295ea_info0.8 \
# --keep /ichrogene/project/temp/gylee/0.GWAS/05.Imputation/2.minimac4/minimac4.1.0/${Sample}/${Sample}.IID \
# --make-pgen \
# --memory ${MEM} \
# --out ${Output_dir}/imputation

#${p2_Dir}/plink2 --pfile ${Output_dir}/imputation \
# --geno 0.02 \
# --make-pgen \
# --memory ${MEM} \
# --update-sex /ichrogene/project/temp/gylee/0.GWAS/KoGES/KNIH_72295ea_info0.8.for_sex.txt \
# --out ${Output_dir}/imputation_g

#${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g \
# --mind 0.02 \
# --make-pgen \
# --memory ${MEM} \
# --out ${Output_dir}/imputation_g_m

#${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m \
# --maf 0.01 \
# --make-pgen \
# --memory ${MEM} \
# --out ${Output_dir}/imputation_g_m_maf

#${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf \
# --hwe 1e-50 \
# --memory ${MEM} \
# --make-pgen \
# --out ${Output_dir}/imputation_g_m_maf_hwe

#${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe \
# --exclude ${Inversion} \
# --indep-pairwise 50 5 0.02 \
# --threads ${Thread} \
# --memory ${MEM} \
# --out ${Output_dir}/${Sample}_indepSNP

#${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe \
# --double-id \
# --extract ${Output_dir}/${Sample}_indepSNP.prune.in \
# --pca approx 10 \
# --memory ${MEM} \
# --set-missing-var-ids @:# \
# --threads ${Thread} \
# --out ${Output_dir}/${Sample}_PCA

#singularity exec --bind /ichrogene/:/ichrogene/ /ichrogene/project/temp/gylee/Singularity/GWAS.sif Rscript /ichrogene/project/temp/gylee/0.GWAS/Code/R/covariate_KOR.R

${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe \
 --glm sex hide-covar \
 --memory ${MEM} \
 --covar ${Output_dir}/covariate_forplink.txt \
 --covar-name PC1-PC10 AGE \
 --pheno ${Output_dir}/${Sample}_pheno.txt \
 --covar-variance-standardize \
 --out ${Output_dir}/${Sample}

cp ${Output_dir}/${Sample}.Pheno.glm.linear ${Output_dir}/${Sample}.glm.repli
#sed -i 's/#CHROM/CHROM/g' ${Output_dir}/${Sample}.Pheno.glm.linear
#sed -i 's/\//:/g' ${Output_dir}/${Sample}.Pheno.glm.linear
#sed -i 's/_/:/g' ${Output_dir}/${Sample}.Pheno.glm.linear

grep -v "^#" ${Output_dir}/imputation_g_m_maf_hwe.pvar | awk '{print $3}' > ${Output_dir}/variant.info
#Bonferroni correction & Make FDR file
singularity exec --bind /ichrogene/:/ichrogene/ /ichrogene/project/temp/gylee/Singularity/GWAS.sif Rscript /ichrogene/project/temp/gylee/0.GWAS/Code/R/Manhattan_plot_continuous.R
