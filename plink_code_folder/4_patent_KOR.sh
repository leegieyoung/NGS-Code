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
PRS_dir="/mnt/nas/gylee/0.GWAS/4.PRS/test"

mkdir -p ${QC_dir}
mkdir -p ${Output_dir} 

export QC_DIR="${QC_dir}/"
export RESULT_DIR="${Output_dir}/"
#echo ${ANA_DIR}
export SAMPLE="${Sample}"

echo ""
echo "----------------------------"
echo " Do        ${CUTOFF}"
echo "----------------------------"
echo ""

rm -rf ${Output_dir}/prune/raw.${CUTOFF}.glm
touch ${Output_dir}/prune/raw.${CUTOFF}.glm

for A in $(cat ${Output_dir}/${CUTOFF}.list)
do
grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid >> ${Output_dir}/prune/raw.${CUTOFF}.glm
done

head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
cat ${Output_dir}/header ${Output_dir}/prune/raw.${CUTOFF}.glm > ${Output_dir}/prune/${CUTOFF}.glm
rm -rf ${Output_dir}/prune/raw.${CUTOFF}.glm


export CUTOFF="${CUTOFF}"
singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker.R
bash /mnt/nas/gylee/0.GWAS/4.PRS/test/calculate_PRS.230713.KOR.sh ${Sample} type2diabet obesity ${CUTOFF}

echo ""
echo "----------------------------"
echo " Finish     ${CUTOFF}"
echo "----------------------------"
echo ""
