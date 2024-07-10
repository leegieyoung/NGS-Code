#!/bin/sh
if [ $# -ne 2 ];then
        echo "Please enter Sample_Name"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
#beagle 진행시Same position 제거를 안했을 경우 사용하는 코드
Sample=$1
Chr=$2
Thread=$3
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

#${p2_Dir}/plink2 --bfile ${QC_dir}/${Sample} \
${p2_Dir}/plink2 --pfile /mnt/nas/gylee/0.GWAS/05.Imputation/2.minimac4/minimac4.1.0/plink/imputation_${Chr} \
 --geno 0.02 \
 --make-pgen \
 --pheno ${QC_dir}/${Sample}_pheno.txt \
 --threads ${Thread} \
 --out ${QC_dir}/${Sample}_g

rm -rf ${QC_dir}/${Sample}_g.pgen
rm -rf ${QC_dir}/${Sample}_g.psam
rm -rf ${QC_dir}/${Sample}_g.pvar

${p2_Dir}/plink2 --bfile ${QC_dir}/${Sample}_g \
 --maf 0.01 \
 --make-psam \
 --threads ${Thread} \
 --out ${QC_dir}/${Sample}_g_maf

rm -rf ${QC_dir}/${Sample}_g.pgen
rm -rf ${QC_dir}/${Sample}_g.psam
rm -rf ${QC_dir}/${Sample}_g.pvar

${p2_Dir}/plink2 --bfile ${QC_dir}/${Sample}_g_maf \
 --hwe 1e-6 \
 --make-psam \
 --threads ${Thread} \
 --out ${QC_dir}/${Sample}_g_maf_hwe

rm -rf ${QC_dir}/${Sample}_g_maf.pgen
rm -rf ${QC_dir}/${Sample}_g_maf.psam
rm -rf ${QC_dir}/${Sample}_g_maf.pvar

echo ""
echo "===================Basic QC End========================"
echo ""

#R pwd
export QC_DIR="${QC_dir}/"
echo ${ANA_DIR}
export SAMPLE="${Sample}"

echo ""
echo "===================Check heterozygosity========================"
echo ""

${p2_Dir}/plink2 --bfile ${QC_dir}/${Sample}_g_maf_hwe \
 --exclude ${Inversion} \
 --indep-pairwise 50 5 0.2 \
 --out ${QC_dir}/${Sample}_indepSNP

${p2_Dir}/plink2 --bfile ${QC_dir}/${Sample}_g_maf_hwe \
 --missing \
 --out ${QC_dir}/Before_miss

${p2_Dir}/plink2 --bfile ${QC_dir}/${Sample}_g_maf_hwe \
 --extract ${QC_dir}/${Sample}_indepSNP.prune.in --het \
 --out ${QC_dir}/Before_het

awk '{print $6"\tPheno"}' ${QC_dir}/Before_miss.smiss > ${QC_dir}/B_F_MISS
paste -d '\t' ${QC_dir}/Before_het.het ${QC_dir}/B_F_MISS > ${QC_dir}/R_BeforeQC_check.txt

singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/het-missing.R

