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
echo ${ANA_DIR}
export SAMPLE="${Sample}"

#cp /mnt/nas/gylee/0.GWAS/05.Imputation/2.minimac4/minimac4.1.0/${Sample}/${Sample}_pheno.txt ${Output_dir}/${Sample}_pheno.txt
#
#${p2_Dir}/plink2 --pfile /mnt/nas/gylee/0.GWAS/1.Input/KNIH_72295ea_info0.8/1.QC/QC_Imputed_KNIH_72295ea_info0.8 \
# --keep /mnt/nas/gylee/0.GWAS/05.Imputation/2.minimac4/minimac4.1.0/${Sample}/${Sample}.IID \
# --make-pgen \
# --memory ${MEM} \
# --out ${Output_dir}/imputation
#
#${p2_Dir}/plink2 --pfile ${Output_dir}/imputation \
# --geno 0.02 \
# --make-pgen \
# --memory ${MEM} \
# --update-sex /mnt/nas/gylee/0.GWAS/1.Input/KNIH_72295ea_info0.8/1.QC/KNIH_72295ea_info0.8.for_sex.txt \
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
#singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/covariate_forSAIGE_KOR.R
#mkdir -p ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LOG
#cp ${Output_dir}/covariate.txt ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/${Sample}.covariate.txt
#
#Log="${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LOG/"

##Divide hwe file to chr
#nohup bash /mnt/nas/gylee/0.GWAS/Code/SAIGE_code_folder/2_divide.chr.sh ${Sample} ${Thread} ${MEM} > ${Log}/plink.log 2>&1 &

##Divide Prune to chr
${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe \
 --extract ${Output_dir}/${Sample}_indepSNP.prune.in \
 --memory ${MEM} \
 --make-bed \
 --out ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC_${Sample}

for A in $(seq 1 22)
do
${p2_Dir}/plink2 --bfile ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC_${Sample} \
 --extract ${Output_dir}/${Sample}_indepSNP.prune.in \
 --memory ${MEM} \
 --make-bed \
 --chr ${A} \
 --out ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${A}.LD_pruned_QC_${Sample}
done

singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/SAIGE/SAIGE_docker_R36.sif bash /mnt/nas/gylee/0.GWAS/Code/SAIGE_code_folder/Step0.createSparseGRM.sh ${Sample} ${Thread}
singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/SAIGE/SAIGE_docker_R36.sif bash /mnt/nas/gylee/0.GWAS/Code/SAIGE_code_folder/Step1.fitNULLGLMM.sh ${Sample} ${Thread}
singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/SAIGE/SAIGE_docker_R36.sif bash /mnt/nas/gylee/0.GWAS/Code/SAIGE_code_folder/Step2.SPAtests_exceptchr2.sh ${Sample} ${Thread}
singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/SAIGE/SAIGE_docker_R36.sif bash /mnt/nas/gylee/0.GWAS/Code/SAIGE_code_folder/Step2.SPAtests_chr2.sh ${Sample} ${Thread}

