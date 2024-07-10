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
CUTOFF=$4
Dir="/mnt/nas/gylee/0.GWAS"
QC_dir="${Dir}/1.Input/${Sample}/1.QC"
MaMi_dir="${Dir}/1.Input/${Sample}/2.MaMi"
p2_Dir="/mnt/nas/gylee/Singurality/plink2"
Inversion="/mnt/nas/gylee/0.GWAS/REFERENCE/inversion.txt"
Output_dir="${Dir}/2.plink_result/${Sample}"
SAIGE_dir="${Dir}/2.SAIGE_result/${Sample}/result"
mkdir -p ${QC_dir}
mkdir -p ${Output_dir} 

export QC_DIR="${QC_dir}/"
export RESULT_DIR="${Output_dir}/"
echo ${ANA_DIR}
export SAMPLE="${Sample}"
export SAIGE_DIR="${SAIGE_dir}/"

##mkdir ${Output_dir}/prune
##${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe \
## --extract ${Output_dir}/prune/${CUTOFF}.list \
## --snps-only just-acgt \
## --make-pgen \
## --out ${Output_dir}/prune/${CUTOFF}.just-acgt
##
##${p2_Dir}/plink2 --pfile ${Output_dir}/prune/${CUTOFF}.just-acgt \
## --indep-pairwise 250 2 0.01 \
## --pheno ${Output_dir}/${Sample}_pheno.txt \
## --memory ${MEM} \
## --out ${Output_dir}/prune/${CUTOFF}
##
rm -rf ${Output_dir}/prune/raw.${CUTOFF}.glm
sed -i -e 's/_/:/g' -e 's/\//:/g' ${Output_dir}/prune/${CUTOFF}.prune.in
rm -rf ${Output_dir}/prune/${CUTOFF}.ori.prune.in
rm -rf ${Output_dir}/prune/raw.${CUTOFF}.ori.prune.in

for A in $(cat ${Output_dir}/prune/${CUTOFF}.prune.in)
do 
grep "${A}" /mnt/nas/gylee/0.GWAS/2.plink_result/Ukb_Train_cliT2D/imputation_g_m_maf_hwe.pvar >> ${Output_dir}/prune/raw.${CUTOFF}.ori.prune.in
done
awk '{print $3}' ${Output_dir}/prune/raw.${CUTOFF}.ori.prune.in > ${Output_dir}/prune/${CUTOFF}.ori.prune.in
rm -rf ${Output_dir}/prune/raw.${CUTOFF}.ori.prune.in
for A in $(cat ${Output_dir}/prune/${CUTOFF}.ori.prune.in)
do
#grep -w "${A}" ${SAIGE_dir}/${Sample}.result >> ${Output_dir}/prune/raw.${CUTOFF}.glm #SAIGE
grep -w "${A}" /mnt/nas/gylee/0.GWAS/2.plink_result/Ukb_Train_cliT2D/firth/${Sample}.Pheno.glm.firth >> ${Output_dir}/prune/raw.${CUTOFF}.ori.glm #plink, firth 
done

#head -n 1 ${SAIGE_dir}/${Sample}.result > ${Output_dir}/header #SAIGE
head -n 1 /mnt/nas/gylee/0.GWAS/2.plink_result/Ukb_Train_cliT2D/firth/${Sample}.Pheno.glm.firth > ${Output_dir}/header #plink, firth

cat ${Output_dir}/header ${Output_dir}/prune/raw.${CUTOFF}.ori.glm > ${Output_dir}/prune/${CUTOFF}.ori.glm
rm -rf ${Output_dir}/prune/raw.${CUTOFF}.ori.glm

rm -rf ${Output_dir}/prune/${CUTOFF}.ori.just-acgt*

export CUTOFF=${CUTOFF}.ori
#singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker_forSAIGE.R
singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker.R #plink, firth

##bash /mnt/nas/gylee/0.GWAS/4.PRS/KoGES/calculate_PRS.230802.KOR_forVali_cliT2D.sh ${Sample} T2D Obesity ${CUTOFF}.ori
##bash /mnt/nas/gylee/0.GWAS/4.PRS/KoGES/calculate_PRS.230811.KOR_forTrain_cliT2D.MH.sh ${Sample} T2D Obesity ${CUTOFF}.ori
bash /mnt/nas/gylee/0.GWAS/4.PRS/test_91584/calculate_PRS_91584.230905_MH.sh ${Sample} type2diabet Obesity ${CUTOFF}
