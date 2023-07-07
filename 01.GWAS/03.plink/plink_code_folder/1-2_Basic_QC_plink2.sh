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

#R pwd
export QC_DIR="${QC_dir}/"
echo ${ANA_DIR}
export SAMPLE="${Sample}"

#1.SampleQC(geno mind impute-sex hwe)
#${p2_Dir}/plink2 --pfile ${QC_dir}/${Sample}_g_m_maf_hwe_het \
# --double-id \
# --pca 20 \
# --set-missing-var-ids @:# \
# --threads ${Thread} \
# --out ${QC_dir}/${Sample}_g_m_maf_hwe_het_PCA
#
#awk '{print $3, $4}' ${QC_dir}/${Sample}_g_m_maf_hwe_het_PCA.eigenvec > ${QC_dir}/PC12
#awk '{print $6}' ${QC_dir}/${Sample}_g_m_maf_hwe_het.fam | sed '1iPheno' > ${QC_dir}/pheno
#sed -i 's/2/Case/g' ${QC_dir}/pheno
#sed -i 's/1/Control/g' ${QC_dir}/pheno
#paste -d ' ' ${QC_dir}/PC12 ${QC_dir}/pheno >  ${QC_dir}/PCA.txt
#singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/PCA.R

echo ""
echo "===================Check IBD========================"
echo ""


#${p2_Dir}/plink2 --pfile ${QC_dir}/${Sample}_g_m_maf_hwe_het \
# --king-cutoff 0.177 \
# --make-bed \
# --pheno ${QC_dir}/${Sample}_pheno.txt \
# --out ${QC_dir}/${Sample}_g_m_maf_hwe_het_KC1st


${p2_Dir}/plink2 --pfile ${QC_dir}/${Sample}_g_m_maf_hwe_het \
 --king-cutoff 0.177 \
 --make-pgen \
 --update-sex /mnt/nas/gylee/0.GWAS/1.Input/Ukb_obesity/1.raw/subset_of_ukb45411_for_sex.txt \
 --out ${QC_dir}/${Sample}_g_m_maf_hwe_het_KC1st

#sh /mnt/nas/gylee/0.GWAS/Code/plink_code_folder/2_ANA_plink2.sh ${Sample} ${Thread}
