#!/bin/sh
if [ $# -ne 5 ];then
        echo "Please enter Sample_Name"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
#beagle 진행시Same position 제거를 안했을 경우 사용하는 코드
Sample=$1
Thread=$2
MEM=$3
CUTOFF=$4
TRAIT=$5

Dir="/ichrogene/project/temp/gylee/0.GWAS"
QC_dir="${Dir}/1.Input/${Sample}/1.QC"
MaMi_dir="${Dir}/1.Input/${Sample}/2.MaMi"
p2_Dir="/ichrogene/project/temp/gylee/Singularity/plink2"
Inversion="/ichrogene/project/temp/gylee/0.GWAS/REFERENCE/inversion.txt"
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
for A in $(cat ${Output_dir}/prune/${CUTOFF}.prune.in)
do
grep -w "${A}" ${SAIGE_dir}/${Sample}.result >> ${Output_dir}/prune/raw.${CUTOFF}.glm #SAIGE
#grep -w "${A}" /ichrogene/project/temp/gylee/0.GWAS/2.plink_result/Ukb_Train_cliT2D/firth/${Sample}.Pheno.glm.firth >> ${Output_dir}/prune/raw.${CUTOFF}.glm #plink, firth 
done

head -n 1 ${SAIGE_dir}/${Sample}.result > ${Output_dir}/header #SAIGE
#head -n 1 /ichrogene/project/temp/gylee/0.GWAS/2.plink_result/Ukb_Train_cliT2D/firth/${Sample}.Pheno.glm.firth > ${Output_dir}/header #plink, firth

cat ${Output_dir}/header ${Output_dir}/prune/raw.${CUTOFF}.glm > ${Output_dir}/prune/${CUTOFF}.glm
rm -rf ${Output_dir}/prune/raw.${CUTOFF}.glm

#rm -rf ${Output_dir}/prune/${CUTOFF}.just-acgt*

export CUTOFF=${CUTOFF}
Rscript /ichrogene/project/temp/gylee/0.GWAS/Code/R/make_marker_forSAIGE.R
#singularity exec --bind /mnt/:/mnt/ /ichrogene/project/temp/gylee/Singularity/GWAS.sif Rscript /ichrogene/project/temp/gylee/0.GWAS/Code/R/make_marker.R #plink, firth

#bash /ichrogene/project/temp/gylee/0.GWAS/4.PRS/Cancer/calculate_PRS.231122_15.6.MH_force_samples.sh ${Sample} ${TRAIT} ${CUTOFF}
bash /ichrogene/project/temp/gylee/0.GWAS/4.PRS/Cancer/calculate_PRS.231127_15.6.MH_force_samples.sh ${Sample} ${TRAIT} ${CUTOFF}
