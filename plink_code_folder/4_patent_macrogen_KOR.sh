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
echo ${ANA_DIR}
export SAMPLE="${Sample}"


#macrogen
rm -rf ${Output_dir}/prune/raw.macrogen.glm
for A in $(cat /mnt/nas/gylee/0.GWAS/1.Input/KNIH_72295ea_info0.8/marker/macrogen_10_KoGES_marker.list)
do
grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid >> ${Output_dir}/prune/raw.macrogen.glm
done

head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
cat ${Output_dir}/header ${Output_dir}/prune/raw.macrogen.glm > ${Output_dir}/prune/macrogen.glm
rm -rf ${Output_dir}/prune/raw.macrogen.glm

export CUTOFF="macrogen"
CUTOFF="macrogen"
singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker.R
bash /mnt/nas/gylee/0.GWAS/4.PRS/test/calculate_PRS.230713.KOR.sh ${Sample} type2diabet obesity ${CUTOFF}

#premium_175
rm -rf ${Output_dir}/prune/raw.premium_175.glm
for A in $(cat /mnt/nas/gylee/0.GWAS/1.Input/KNIH_72295ea_info0.8/marker/premium_175_KoGES_marker.list)
do
grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid >> ${Output_dir}/prune/raw.premium_175.glm
done

head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
cat ${Output_dir}/header ${Output_dir}/prune/raw.premium_175.glm > ${Output_dir}/prune/premium_175.glm
rm -rf ${Output_dir}/prune/raw.premium_175.glm

export CUTOFF="premium_175"
CUTOFF="premium_175"
singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker.R
bash /mnt/nas/gylee/0.GWAS/4.PRS/test/calculate_PRS.230713.KOR.sh ${Sample} type2diabet obesity ${CUTOFF}


