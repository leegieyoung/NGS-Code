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
#echo ${ANA_DIR}
export SAMPLE="${Sample}"

#awk '{print $3}' ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid.FDR > ${Output_dir}/${Sample}.snplist
#awk '{print $3}' ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid.e-4 > ${Output_dir}/${Sample}.e-4.snplist
#awk '{print $3}' ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid.e-4 > ${Output_dir}/${Sample}.e-4.snplist

#Premium markers 173 T2D
#awk '{print $1}' ${Dir}/4.PRS/test/type2diabet/rsID2ChrPosiEA_type2diabet_marker.txt > ${Output_dir}/prune/premium_173.T2D.txt
#
###FDR.min0.02 Marker
#echo ""
#echo "----------------------------"
#echo " Do        FDR.min0.02"
#echo "----------------------------"
#echo ""
#cat ${Output_dir}/prune/premium_173.T2D.txt ${Dir}/2.plink_result/Imputed_Ukb_obesity_BMI_WC_premium/prune/FDR.min0.02.prune.in | sort | uniq > ${Output_dir}/prune/Obesity.FDR.min0.02_173.T2D.txt 
#
##FDR.min0.02 Markers
#rm -rf ${Output_dir}/prune/raw.Obesity.FDR.min0.02_173.T2D.glm
#touch ${Output_dir}/prune/raw.Obesity.FDR.min0.02_173.T2D.glm
#
##for A in $(cat ${Output_dir}/prune/Obesity.FDR.min0.02_173.T2D.txt)
##do
##grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid >> ${Output_dir}/prune/raw.Obesity.FDR.min0.02_173.T2D.glm
##done
#
#head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
#cat ${Output_dir}/header ${Output_dir}/prune/raw.Obesity.FDR.min0.02_173.T2D.glm  > ${Output_dir}/prune/Obesity.FDR.min0.02_173.T2D.glm
#rm -rf ${Output_dir}/prune/raw.Obesity.FDR.min0.02_173.T2D.glm
#
##export CUTOFF="FDR.min0.02_173.T2D"
##singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker_Obe_T2D.R
#
##nohup bash ${PRS_dir}/calculate_PRS.FDR.min0.02_173.T2D.sh ${Sample} type2diabet obesity &
#
#echo ""
#echo "----------------------------"
#echo " Finish     FDR.min0.02"
#echo "----------------------------"
#echo ""
#
###r2 < 0.02, FDR Markers
#echo ""
#echo "----------------------------"
#echo " Do        FDR.min0.02"
#echo "----------------------------"
#echo ""
#cat ${Output_dir}/prune/premium_173.T2D.txt ${Dir}/2.plink_result/Imputed_Ukb_obesity_BMI_WC_premium/prune/FDR.min0.02.prune.in | sort | uniq > ${Output_dir}/prune/Obesity.FDR.min0.02.txt 
#
###FDR.min0.02 Markers
#rm -rf ${Output_dir}/prune/raw.Obesity.FDR.min0.02.glm
#touch ${Output_dir}/prune/raw.Obesity.FDR.min0.02.glm
#
##for A in $(cat ${Output_dir}/prune/Obesity.FDR.min0.02.txt)
##do
##grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid >> ${Output_dir}/prune/raw.Obesity.FDR.min0.02.glm
##done
#
#head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
#cat ${Output_dir}/header ${Output_dir}/prune/raw.Obesity.FDR.min0.02.glm  > ${Output_dir}/prune/Obesity.FDR.min0.02.glm
#rm -rf ${Output_dir}/prune/raw.Obesity.FDR.min0.02.glm
#
##export CUTOFF="FDR.min0.02"
##singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker_Obe_T2D.R
#
##nohup bash ${PRS_dir}/calculate_PRS.FDR.min0.02.sh ${Sample} type2diabet obesity &
#
#echo ""
#echo "----------------------------"
#echo " Finish     FDR.min0.02"
#echo "----------------------------"
#echo ""
#
#
##e-4.min0.02 Marker
#echo ""
#echo "----------------------------"
#echo " Do        e-4.min0.02"
#echo "----------------------------"
#echo ""
#cat ${Output_dir}/prune/premium_173.T2D.txt ${Dir}/2.plink_result/Imputed_Ukb_obesity_BMI_WC_premium/prune/e-4.min0.02.prune.in | sort | uniq > ${Output_dir}/prune/Obesity.e-4.min0.02_173.T2D.txt
#rm -rf ${Output_dir}/prune/raw.Obesity.e-4.min0.02_173.T2D.glm
#touch ${Output_dir}/prune/raw.Obesity.e-4.min0.02_173.T2D.glm
#
#for A in $(cat ${Output_dir}/prune/Obesity.e-4.min0.02_173.T2D.txt)
#do
#grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid >> ${Output_dir}/prune/raw.Obesity.e-4.min0.02_173.T2D.glm
#done
#
#head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
#cat ${Output_dir}/header ${Output_dir}/prune/raw.Obesity.e-4.min0.02_173.T2D.glm > ${Output_dir}/prune/Obesity.e-4.min0.02_173.T2D.glm
#rm -rf ${Output_dir}/prune/raw.Obesity.e-4.min0.02_173.T2D.glm
#
#export CUTOFF="e-4.min0.02_173.T2D"
#singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker_Obe_T2D.R
#
#nohup bash ${PRS_dir}/calculate_PRS.e-4.min0.02_173.T2D.sh ${Sample} type2diabet obesity &
#
#echo ""
#echo "----------------------------"
#echo " Finish     e-4.min0.02"
#echo "----------------------------"
#echo ""
#
##e-4.min0.02 Marker
echo ""
echo "----------------------------"
echo " Do        e-4.min0.02"
echo "----------------------------"
echo ""

rm -rf ${Output_dir}/prune/raw.Obesity.e-4.min0.02.glm
touch ${Output_dir}/prune/raw.Obesity.e-4.min0.02.glm

for A in $(cat ${Output_dir}/prune/e-4.min0.02.prune.in)
do
grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid >> ${Output_dir}/prune/raw.Obesity.e-4.min0.02.glm
done

head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
cat ${Output_dir}/header ${Output_dir}/prune/raw.Obesity.e-4.min0.02.glm > ${Output_dir}/prune/Obesity.e-4.min0.02.glm
rm -rf ${Output_dir}/prune/raw.Obesity.e-4.min0.02.glm

export CUTOFF="e-4.min0.02"
singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker_Obe_T2D.R

nohup bash ${PRS_dir}/calculate_PRS.e-4.min0.02.sh ${Sample} type2diabet obesity &

echo ""
echo "----------------------------"
echo " Finish     e-4.min0.02"
echo "----------------------------"
echo ""


##e-4.min0.02 Marker
#echo ""
#echo "----------------------------"
#echo " Do        e-4.min0.02"
#echo "----------------------------"
#echo ""
#cat ${Output_dir}/prune/premium_174.T2D.txt ${Dir}/2.plink_result/Imputed_Ukb_obesity_BMI_WC_premium/prune/e-4.min0.02.prune.in | sort | uniq > ${Output_dir}/prune/Obesity.e-4.min0.02_173.T2D.txt
#rm -rf ${Output_dir}/prune/raw.Obesity.2e-4.min0.02_173.T2D.glm
#touch ${Output_dir}/prune/raw.Obesity.e-4.min0.02_173.T2D.glm
#
##for A in $(cat ${Output_dir}/prune/Obesity.e-4.min0.02_173.T2D.txt)
##do
##grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid >> ${Output_dir}/prune/raw.Obesity.e-4.min0.02_173.T2D.glm
##done
#
#head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
#cat ${Output_dir}/header ${Output_dir}/prune/raw.Obesity.e-4.min0.02_173.T2D.glm > ${Output_dir}/prune/Obesity.e-4.min0.02_173.T2D.glm
#rm -rf ${Output_dir}/prune/raw.Obesity.e-4.min0.02_173.T2D.glm
#
#export CUTOFF="e-4.min0.02_173.T2D"
#singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker_Obe_T2D.R
#
#echo ""
#echo "----------------------------"
#echo " Finish     e-4.min0.02"
#echo "----------------------------"
#echo ""
#
##e-4.min0.02 Marker
#echo ""
#echo "----------------------------"
#echo " Do        e-4.min0.02"
#echo "----------------------------"
#echo ""
#cat ${Output_dir}/prune/premium_173.T2D.txt ${Dir}/2.plink_result/Imputed_Ukb_obesity_BMI_WC_premium/prune/e-4.min0.02.prune.in | sort | uniq > ${Output_dir}/prune/Obesity.e-4.min0.02.txt
#rm -rf ${Output_dir}/prune/raw.Obesity.e-4.min0.02.glm
#touch ${Output_dir}/prune/raw.Obesity.e-4.min0.02.glm
#
##for A in $(cat ${Output_dir}/prune/Obesity.e-4.min0.02.txt)
##do
##grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid >> ${Output_dir}/prune/raw.Obesity.e-4.min0.02.glm
##done
#
#head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
#cat ${Output_dir}/header ${Output_dir}/prune/raw.Obesity.e-4.min0.02.glm > ${Output_dir}/prune/Obesity.e-4.min0.02.glm
#rm -rf ${Output_dir}/prune/raw.Obesity.e-4.min0.02.glm
#
#export CUTOFF="e-4.min0.02"
#singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker_Obe_T2D.R
#
#echo ""
#echo "----------------------------"
#echo " Finish     e-4.min0.02"
#echo "----------------------------"
#echo ""
