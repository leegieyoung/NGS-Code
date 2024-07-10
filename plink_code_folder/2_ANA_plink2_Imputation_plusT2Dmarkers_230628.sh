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


#Premium markers 173 T2D
awk '{print $1}' /mnt/nas/gylee/0.GWAS/4.PRS/test/type2diabet/rsID2ChrPosiEA_type2diabet_marker.txt > ${Output_dir}/prune/premium_173_T2D.txt

#Marker
rm -rf ${Output_dir}/prune/premium_173_T2D.glm
for A in $(cat ${Output_dir}/prune/premium_173_T2D.txt)
do
grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid >> ${Output_dir}/prune/premium_173_T2D.glm
done
echo "-------"
echo "${Sample} in 173T2D markers"
wc -l ${Output_dir}/prune/premium_173_T2D.glm
echo "-------"

mkdir ${Output_dir}/prune
${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe \
 --extract ${Output_dir}/${Sample}.snplist \
 --snps-only just-acgt \
 --make-pgen \
 --out ${Output_dir}/prune/just-acgt
 
${p2_Dir}/plink2 --pfile ${Output_dir}/prune/just-acgt \
 --indep-pairwise 1000kb 0.2 \
 --pheno ${Output_dir}/${Sample}_pheno.txt \
 --memory ${MEM} \
 --out ${Output_dir}/prune/min0.2

rm -rf ${Output_dir}/prune/raw_prune.glm
touch ${Output_dir}/prune/raw_prune.glm
 
for A in $(cat ${Output_dir}/prune/min0.2.prune.in)
do
grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid.FDR >> ${Output_dir}/prune/raw_prune.glm
done
head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
cat ${Output_dir}/header ${Output_dir}/prune/premium_173_T2D.glm ${Output_dir}/prune/raw_prune.glm | sort | uniq > ${Output_dir}/prune/prune.glm
 
${p2_Dir}/plink2 --pfile ${Output_dir}/prune/just-acgt \
 --indep-pairwise 1000kb 0.02 \
 --pheno ${Output_dir}/${Sample}_pheno.txt \
 --memory ${MEM} \
 --out ${Output_dir}/prune/min0.02

rm -rf ${Output_dir}/prune/raw_min0.02.prune.glm
touch ${Output_dir}/prune/raw_min0.02prune.glm
 
for A in $(cat ${Output_dir}/prune/min0.02.prune.in)
do
grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid.FDR >> ${Output_dir}/prune/raw_min0.02.prune.glm
done
head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
cat ${Output_dir}/header ${Output_dir}/prune/premium_173_T2D.glm ${Output_dir}/prune/raw_min0.02.prune.glm | sort | uniq > ${Output_dir}/prune/min0.02.prune.glm
#singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/3_VEP_uniq_gene.R

#Suggest
mkdir ${Output_dir}/prune
${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe \
 --extract ${Output_dir}/${Sample}.suggest.snplist \
 --snps-only just-acgt \
 --make-pgen \
 --out ${Output_dir}/prune/suggest.just-acgt

${p2_Dir}/plink2 --pfile ${Output_dir}/prune/suggest.just-acgt \
 --indep-pairwise 1000kb 0.2 \
 --pheno ${Output_dir}/${Sample}_pheno.txt \
 --memory ${MEM} \
 --out ${Output_dir}/prune/suggest.min0.2

rm -rf ${Output_dir}/prune/raw_prune.suggest.glm
touch ${Output_dir}/prune/raw_prune.suggest.glm

for A in $(cat ${Output_dir}/prune/suggest.min0.2.prune.in)
do
grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid.suggest >> ${Output_dir}/prune/raw_prune.suggest.glm
done
head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
cat ${Output_dir}/header ${Output_dir}/prune/premium_173_T2D.glm ${Output_dir}/prune/raw_prune.suggest.glm | sort | uniq  > ${Output_dir}/prune/prune.suggest.glm

#1e-4
mkdir ${Output_dir}/prune
${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe \
 --extract ${Output_dir}/${Sample}.1e-4.snplist \
 --snps-only just-acgt \
 --make-pgen \
 --out ${Output_dir}/prune/1e-4.just-acgt

${p2_Dir}/plink2 --pfile ${Output_dir}/prune/1e-4.just-acgt \
 --indep-pairwise 1000kb 0.2 \
 --pheno ${Output_dir}/${Sample}_pheno.txt \
 --memory ${MEM} \
 --out ${Output_dir}/prune/1e-4.min0.2

rm -rf ${Output_dir}/prune/raw_prune.1e-4.glm
touch ${Output_dir}/prune/raw_prune.1e-4.glm

for A in $(cat ${Output_dir}/prune/1e-4.min0.2.prune.in)
do
grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid.1e-4 >> ${Output_dir}/prune/raw_prune.1e-4.glm
done
head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
cat ${Output_dir}/header ${Output_dir}/prune/premium_173_T2D.glm ${Output_dir}/prune/raw_prune.1e-4.glm | sort | uniq > ${Output_dir}/prune/prune.1e-4.glm

singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker.R

