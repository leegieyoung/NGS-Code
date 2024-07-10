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

#1.SampleQC(geno mind impute-sex hwe)
mkdir -p ${QC_dir}
mkdir -p ${MaMi_dir}

echo ""
echo "===================PWD & Sample names========================"
echo ""

echo ${QC_dir}
echo ${Sample}
echo ""
echo "===================PWD & Sample names========================"
echo ""

echo ""
echo "===================Basic QC Start========================"
echo ""

${p2_Dir}/plink2 --bfile ${QC_dir}/${Sample} \
 --geno 0.02 \
 --make-bed \
 --pheno ${QC_dir}/${Sample}_pheno.txt \
 --threads ${Thread} \
 --out ${QC_dir}/${Sample}_g


${p2_Dir}/plink2 --bfile ${QC_dir}/${Sample}_g \
 --mind 0.02 \
 --make-bed \
 --threads ${Thread} \
 --out ${QC_dir}/${Sample}_g_m

${p2_Dir}/plink2 --bfile ${QC_dir}/${Sample}_g_m \
 --maf 0.01 \
 --make-bed \
 --threads ${Thread} \
 --out ${QC_dir}/${Sample}_g_m_maf

${p2_Dir}/plink2 --bfile ${QC_dir}/${Sample}_g_m_maf \
 --hwe 1e-6 \
 --make-bed \
 --threads ${Thread} \
 --out ${QC_dir}/${Sample}_g_m_maf_hwe

echo ""
echo "===================Basic QC End========================"
echo ""

echo ""
echo "===================Check MAF========================"
echo ""

${p2_Dir}/plink2 --bfile ${QC_dir}/${Sample}_g_m_maf_hwe \
 --freq \
 --out ${QC_dir}/${Sample}_g_m_maf_hwe_freq

#R pwd
export QC_DIR="${QC_dir}/"
echo ${ANA_DIR}
export SAMPLE="${Sample}"

singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/MAF_check.R

echo ""
echo "===================Check heterozygosity========================"
echo ""

${p2_Dir}/plink2 --bfile ${QC_dir}/${Sample}_g_m_maf_hwe \
 --exclude ${Inversion} \
 --indep-pairwise 50 5 0.2 \
 --out ${QC_dir}/${Sample}_indepSNP


${p2_Dir}/plink2 --bfile ${QC_dir}/${Sample}_g_m_maf_hwe \
 --extract ${QC_dir}/${Sample}_indepSNP.prune.in --het --missing \
 --out ${QC_dir}/${Sample}_R_check

awk '{print $5"\tPheno"}' ${QC_dir}/${Sample}_R_check.smiss > ${QC_dir}/F_MISS
paste -d '\t' ${QC_dir}/${Sample}_R_check.het ${QC_dir}/F_MISS > ${QC_dir}/R_check.txt 

singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/heterozygosity_outliers_list.R

sed 's/"// g' ${QC_dir}/${Sample}_fail-het-qc.txt | awk '{print$1, $2}'> ${QC_dir}/${Sample}_het_fail_ind.txt

${p2_Dir}/plink2 --bfile ${QC_dir}/${Sample}_g_m_maf_hwe \
 --remove ${QC_dir}/${Sample}_het_fail_ind.txt \
 --mind 0.02 \
 n--make-bed \
 --out ${QC_dir}/${Sample}_g_m_maf_hwe_het

echo ""
echo "===================Check PCA========================"
echo ""

${p2_Dir}/plink2 --bfile ${QC_dir}/${Sample}_g_m_maf_hwe_het \
 --double-id \
 --pca 10 \
 --set-missing-var-ids @:# \
 --threads ${Thread} \
 --out ${QC_dir}/${Sample}_g_m_maf_hwe_het_PCA

awk '{print $3, $4}' ${QC_dir}/${Sample}_g_m_maf_hwe_het_PCA.eigenvec > ${QC_dir}/PC12
awk '{print $6}' ${QC_dir}/${Sample}_g_m_maf_hwe_het.fam | sed '1iPheno' > ${QC_dir}/pheno
sed -i 's/2/Case/g' ${QC_dir}/pheno
sed -i 's/1/Control/g' ${QC_dir}/pheno
paste -d ' ' ${QC_dir}/PC12 ${QC_dir}/pheno >  ${QC_dir}/PCA.txt
singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/PCA.R

#QC check

${p2_Dir}/plink2 --bfile ${QC_dir}/${Sample} \
 --missing \
 --out ${QC_dir}/Before_miss

${p2_Dir}/plink2 --bfile ${QC_dir}/${Sample} \
 --extract ${QC_dir}/${Sample}_indepSNP.prune.in --het \
 --out ${QC_dir}/Before_het

awk '{print $6"\tPheno"}' ${QC_dir}/Before_miss.smiss > ${QC_dir}/B_F_MISS
paste -d '\t' ${QC_dir}/Before_het.het ${QC_dir}/B_F_MISS > ${QC_dir}/R_BeforeQC_check.txt


${p2_Dir}/plink2 --bfile ${QC_dir}/${Sample}_g_m_maf_hwe_het \
 --missing \
 --out ${QC_dir}/After_miss

${p2_Dir}/plink2 --bfile ${QC_dir}/${Sample}_g_m_maf_hwe_het \
 --extract ${QC_dir}/${Sample}_indepSNP.prune.in --het \
 --out ${QC_dir}/After_het

awk '{print $6"\tPheno"}' ${QC_dir}/After_miss.smiss > ${QC_dir}/A_F_MISS
paste -d '\t' ${QC_dir}/After_het.het ${QC_dir}/A_F_MISS > ${QC_dir}/R_AfterQC_check.txt

singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/het-missing.R

echo ""
echo "===================Check IBD========================"
echo ""


#${p2_Dir}/plink2 --bfile ${QC_dir}/${Sample}_g_m_maf_hwe_het \
# --king-cutoff 0.177 \
# --make-bed \
# --pheno ${QC_dir}/${Sample}_pheno.txt \
# --out ${QC_dir}/${Sample}_g_m_maf_hwe_het_KC1st


${p2_Dir}/plink2 --bfile ${QC_dir}/${Sample}_g_m_maf_hwe_het \
 --king-cutoff 0.177 \
 --make-bed \
 --pheno ${QC_dir}/${Sample}_pheno.txt \
 --update-sex /mnt/nas/gylee/0.GWAS/1.Input/Ukb_obesity/1.raw/subset_of_ukb45411_for_sex.txt \
 --out ${QC_dir}/${Sample}_g_m_maf_hwe_het_KC1st

sh /mnt/nas/gylee/0.GWAS/Code/plink_code_folder/2_ANA_plink2.sh ${Sample} ${Thread}
