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
PRS_dir="/mnt/nas/gylee/0.GWAS/4.PRS/test"

mkdir -p ${QC_dir}
mkdir -p ${Output_dir} 

export QC_DIR="${QC_dir}/"
export RESULT_DIR="${Output_dir}/"
export SAMPLE="${Sample}"

echo ""
echo "----------------------------"
echo " Do        UKB_T2D_suggest.min0.02_Obe_suggest.min0.02"
echo "----------------------------"
echo ""
cat ${Output_dir}/prune/suggest.min0.02.prune.in ${Dir}/2.plink_result/Imputed_Ukb_obesity_BMI_WC_premium/prune/suggest.min0.02.prune.in | sort | uniq > ${Output_dir}/prune/UKB_T2D_suggest.min0.02_Obe_suggest.min0.02.txt
rm -rf ${Output_dir}/prune/raw.UKB_T2D_suggest.min0.02_Obe_suggest.min0.02.glm
touch ${Output_dir}/prune/raw.UKB_T2D_suggest.min0.02_Obe_suggest.min0.02.glm

for A in $(cat ${Output_dir}/prune/UKB_T2D_suggest.min0.02_Obe_suggest.min0.02.txt)
do
grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid >> ${Output_dir}/prune/raw.UKB_T2D_suggest.min0.02_Obe_suggest.min0.02.glm
done

head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
cat ${Output_dir}/header ${Output_dir}/prune/raw.UKB_T2D_suggest.min0.02_Obe_suggest.min0.02.glm > ${Output_dir}/prune/UKB_T2D_suggest.min0.02_Obe_suggest.min0.02.glm
rm -rf ${Output_dir}/prune/raw.UKB_T2D_suggest.min0.02_Obe_suggest.min0.02.glm

export CUTOFF="UKB_T2D_suggest.min0.02_Obe_suggest.min0.02"

singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker.R
bash ${PRS_dir}/calculate_PRS.230707.sh ${Sample} type2diabet obesity ${CUTOFF} 

echo ""
echo "----------------------------"
echo " Finish     UKB_T2D_suggest.min0.02_Obe_suggest.min0.02"
echo "----------------------------"
echo ""

