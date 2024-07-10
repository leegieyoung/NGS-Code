#!/bin/sh
if [ $# -ne 6 ];then
        echo "Please enter Sample_Name"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
#beagle 진행시Same position 제거를 안했을 경우 사용하는 코드
Sample=$1
Thread=$2
MEM=$3
CUTOFF=$4
Trait1=$5
TrainVali=$6
Dir="/ichrogene/project/temp/gylee/0.GWAS"
QC_dir="${Dir}/1.Input/${Sample}/1.QC"
MaMi_dir="${Dir}/1.Input/${Sample}/2.MaMi"
p2_Dir="/ichrogene/project/temp/gylee/Singurality/plink2"
Inversion="/ichrogene/project/temp/gylee/0.GWAS/REFERENCE/inversion.txt"
Output_dir="${Dir}/2.plink_result/${TrainVali}/${Sample}"
SAIGE_dir="${Dir}/2.SAIGE_result/${Sample}/result"

mkdir -p ${QC_dir}
mkdir -p ${Output_dir}

export QC_DIR="${QC_dir}/"
export RESULT_DIR="${Output_dir}/"
echo ${ANA_DIR}
export SAMPLE="${Sample}"
export plink_DIR="${plink_dir}/"

#mkdir ${Output_dir}/prune
#${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe \
# --extract ${Output_dir}/prune/${CUTOFF}.list \
# --snps-only just-acgt \
# --make-pgen \
# --out ${Output_dir}/prune/${CUTOFF}.just-acgt
#
#${p2_Dir}/plink2 --pfile ${Output_dir}/prune/${CUTOFF}.just-acgt \
# --indep-pairwise 250 2 0.01 \
# --pheno ${Output_dir}/${Sample}_pheno.txt \
# --memory ${MEM} \
# --out ${Output_dir}/prune/${CUTOFF}
#
rm -rf ${Output_dir}/prune/raw.${CUTOFF}.linear
for A in $(cat ${Output_dir}/prune/${CUTOFF}.prune.in)
do
grep -w "${A}" ${Output_dir}/${Sample}.Pheno.glm.linear >> ${Output_dir}/prune/raw.${CUTOFF}.linear
done

head -n 2 ${Output_dir}/${Sample}.Pheno.glm.linear > ${Output_dir}/header
cat ${Output_dir}/header ${Output_dir}/prune/raw.${CUTOFF}.linear > ${Output_dir}/prune/${CUTOFF}.linear
rm -rf ${Output_dir}/prune/raw.${CUTOFF}.linear

rm -rf ${Output_dir}/prune/${CUTOFF}.just-acgt*

export CUTOFF=${CUTOFF}
#singularity exec --bind /mnt/:/mnt/ /ichrogene/project/temp/gylee/Singurality/GWAS.sif Rscript /ichrogene/project/temp/gylee/0.GWAS/Code/R/make_marker.linear.R
Rscript /ichrogene/project/temp/gylee/0.GWAS/Code/R/make_marker.linear.R
#bash /ichrogene/project/temp/gylee/1.GWAS/4.PRS/KoGES/calculate_PRS.230802.KOR_forVali_cliT2D.sh ${Sample} T2D Obesity ${CUTOFF}
#bash /ichrogene/project/temp/gylee/0.GWAS/4.PRS/KoGES/calculate_PRS.230811.KOR_forTrain_cliT2D.MH.sh ${Sample} T2D Obesity ${CUTOFF}
bash /ichrogene/project/temp/gylee/0.GWAS/4.PRS/KoGES/calculate_PRS.240326.KOR_forTrain_continuous.sh ${Sample} ${Trait1} ${CUTOFF} ${TrainVali}
bash /ichrogene/project/temp/gylee/0.GWAS/4.PRS/KoGES/calculate_PRS.240226.KOR_forVali_continuous.sh ${Sample} ${Trait1} ${CUTOFF} ${TrainVali}

