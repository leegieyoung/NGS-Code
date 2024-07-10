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
#awk '{print $3}' ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid.suggest > ${Output_dir}/${Sample}.suggest.snplist
#awk '{print $3}' ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid.1e-4 > ${Output_dir}/${Sample}.1e-4.snplist

#Premium markers 173 T2D
#awk '{print $1}' ${Dir}/4.PRS/test/type2diabet/rsID2ChrPosiEA_type2diabet_marker.txt > ${Output_dir}/prune/premium_173.T2D.txt
#
###FDR.min0.2 Marker
#echo ""
#echo "----------------------------"
#echo " Do        FDR.min0.2"
#echo "----------------------------"
#echo ""
#cat ${Output_dir}/prune/premium_173.T2D.txt ${Dir}/2.plink_result/Imputed_Ukb_obesity_BMI_WC_premium/prune/FDR.min0.2.prune.in | sort | uniq > ${Output_dir}/prune/Obesity.FDR.min0.2_173.T2D.txt 
#
##FDR.min0.2 Markers
#rm -rf ${Output_dir}/prune/raw.Obesity.FDR.min0.2_173.T2D.glm
#touch ${Output_dir}/prune/raw.Obesity.FDR.min0.2_173.T2D.glm
#
##for A in $(cat ${Output_dir}/prune/Obesity.FDR.min0.2_173.T2D.txt)
##do
##grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid >> ${Output_dir}/prune/raw.Obesity.FDR.min0.2_173.T2D.glm
##done
#
#head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
#cat ${Output_dir}/header ${Output_dir}/prune/raw.Obesity.FDR.min0.2_173.T2D.glm  > ${Output_dir}/prune/Obesity.FDR.min0.2_173.T2D.glm
#rm -rf ${Output_dir}/prune/raw.Obesity.FDR.min0.2_173.T2D.glm
#
##export CUTOFF="FDR.min0.2_173.T2D"
##singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker_Obe_T2D.R
#
##nohup bash ${PRS_dir}/calculate_PRS.FDR.min0.2_173.T2D.sh ${Sample} type2diabet obesity &
#
#echo ""
#echo "----------------------------"
#echo " Finish     FDR.min0.2"
#echo "----------------------------"
#echo ""
#
###r2 < 0.2, FDR Markers
#echo ""
#echo "----------------------------"
#echo " Do        FDR.min0.2"
#echo "----------------------------"
#echo ""
#cat ${Output_dir}/prune/premium_173.T2D.txt ${Dir}/2.plink_result/Imputed_Ukb_obesity_BMI_WC_premium/prune/FDR.min0.2.prune.in | sort | uniq > ${Output_dir}/prune/Obesity.FDR.min0.2_173.T2D.txt 
#
###FDR.min0.2 Markers
#rm -rf ${Output_dir}/prune/raw.Obesity.FDR.min0.2_173.T2D.glm
#touch ${Output_dir}/prune/raw.Obesity.FDR.min0.2_173.T2D.glm
#
##for A in $(cat ${Output_dir}/prune/Obesity.FDR.min0.2_173.T2D.txt)
##do
##grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid >> ${Output_dir}/prune/raw.Obesity.FDR.min0.2_173.T2D.glm
##done
#
#head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
#cat ${Output_dir}/header ${Output_dir}/prune/raw.Obesity.FDR.min0.2_173.T2D.glm  > ${Output_dir}/prune/Obesity.FDR.min0.2_173.T2D.glm
#rm -rf ${Output_dir}/prune/raw.Obesity.FDR.min0.2_173.T2D.glm
#
##export CUTOFF="FDR.min0.2_173.T2D"
##singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker_Obe_T2D.R
#
##nohup bash ${PRS_dir}/calculate_PRS.FDR.min0.2_173.T2D.sh ${Sample} type2diabet obesity &
#
#echo ""
#echo "----------------------------"
#echo " Finish     FDR.min0.2"
#echo "----------------------------"
#echo ""
#
#
##Suggest.min0.2 Marker
#echo ""
#echo "----------------------------"
#echo " Do        Suggest.min0.2"
#echo "----------------------------"
#echo ""
#cat ${Output_dir}/prune/premium_173.T2D.txt ${Dir}/2.plink_result/Imputed_Ukb_obesity_BMI_WC_premium/prune/FDR.min0.2.prune.in | sort | uniq > ${Output_dir}/prune/Obesity.FDR.min0.2_173.T2D.txt
#rm -rf ${Output_dir}/prune/raw.Obesity.FDR.min0.2_173.T2D.glm
#touch ${Output_dir}/prune/raw.Obesity.FDR.min0.2_173.T2D.glm
#
#for A in $(cat ${Output_dir}/prune/Obesity.FDR.min0.2_173.T2D.txt)
#do
#grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid >> ${Output_dir}/prune/raw.Obesity.FDR.min0.2_173.T2D.glm
#done
#
#head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
#cat ${Output_dir}/header ${Output_dir}/prune/raw.Obesity.FDR.min0.2_173.T2D.glm > ${Output_dir}/prune/Obesity.FDR.min0.2_173.T2D.glm
#rm -rf ${Output_dir}/prune/raw.Obesity.FDR.min0.2_173.T2D.glm
#
#export CUTOFF="FDR.min0.2_173.T2D"
#singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker_Obe_T2D.R
#
#nohup bash ${PRS_dir}/calculate_PRS.FDR.min0.2_173.T2D.sh ${Sample} type2diabet obesity &
#
#echo ""
#echo "----------------------------"
#echo " Finish     Suggest.min0.2"
#echo "----------------------------"
#echo ""
#
##Suggest.min0.2 Marker
echo ""
echo "----------------------------"
echo " Do        Suggest.min0.2"
echo "----------------------------"
echo ""
#cat ${Output_dir}/prune/premium_173.T2D.txt ${Dir}/2.plink_result/Imputed_Ukb_obesity_BMI_WC_premium/prune/FDR.min0.2.prune.in | sort | uniq > ${Output_dir}/prune/Obesity.FDR.min0.2_Only_Obesity.txt
rm -rf ${Output_dir}/prune/raw.Obesity.FDR.min0.2_Only_Obesity.glm
touch ${Output_dir}/prune/raw.Obesity.FDR.min0.2_Only_Obesity.glm

for A in $(cat ${Dir}/2.plink_result/Imputed_Ukb_obesity_BMI_WC_premium/prune/FDR.min0.2.prune.in)
do
grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid >> ${Output_dir}/prune/raw.Obesity.FDR.min0.2_Only_Obesity.glm
done

head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
cat ${Output_dir}/header ${Output_dir}/prune/raw.Obesity.FDR.min0.2_Only_Obesity.glm > ${Output_dir}/prune/Obesity.FDR.min0.2_Only_Obesity.glm
rm -rf ${Output_dir}/prune/raw.Obesity.FDR.min0.2_Only_Obesity.glm

export CUTOFF="FDR.min0.2_Only_Obesity"
singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker_Obe_T2D.R

nohup bash ${PRS_dir}/calculate_PRS.FDR.min0.2_Only_Obesity.sh ${Sample} type2diabet obesity &

echo ""
echo "----------------------------"
echo " Finish     Suggest.min0.2"
echo "----------------------------"
echo ""


##1e-4.min0.2 Marker
#echo ""
#echo "----------------------------"
#echo " Do        1e-4.min0.2"
#echo "----------------------------"
#echo ""
#cat ${Output_dir}/prune/premium_174.T2D.txt ${Dir}/2.plink_result/Imputed_Ukb_obesity_BMI_WC_premium/prune/1e-4.min0.2.prune.in | sort | uniq > ${Output_dir}/prune/Obesity.1e-4.min0.2_173.T2D.txt
#rm -rf ${Output_dir}/prune/raw.Obesity.2e-4.min0.2_173.T2D.glm
#touch ${Output_dir}/prune/raw.Obesity.1e-4.min0.2_173.T2D.glm
#
##for A in $(cat ${Output_dir}/prune/Obesity.1e-4.min0.2_173.T2D.txt)
##do
##grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid >> ${Output_dir}/prune/raw.Obesity.1e-4.min0.2_173.T2D.glm
##done
#
#head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
#cat ${Output_dir}/header ${Output_dir}/prune/raw.Obesity.1e-4.min0.2_173.T2D.glm > ${Output_dir}/prune/Obesity.1e-4.min0.2_173.T2D.glm
#rm -rf ${Output_dir}/prune/raw.Obesity.1e-4.min0.2_173.T2D.glm
#
#export CUTOFF="1e-4.min0.2_173.T2D"
#singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker_Obe_T2D.R
#
#echo ""
#echo "----------------------------"
#echo " Finish     1e-4.min0.2"
#echo "----------------------------"
#echo ""
#
##1e-4.min0.2 Marker
#echo ""
#echo "----------------------------"
#echo " Do        1e-4.min0.2"
#echo "----------------------------"
#echo ""
#cat ${Output_dir}/prune/premium_173.T2D.txt ${Dir}/2.plink_result/Imputed_Ukb_obesity_BMI_WC_premium/prune/1e-4.min0.2.prune.in | sort | uniq > ${Output_dir}/prune/Obesity.1e-4.min0.2_173.T2D.txt
#rm -rf ${Output_dir}/prune/raw.Obesity.1e-4.min0.2_173.T2D.glm
#touch ${Output_dir}/prune/raw.Obesity.1e-4.min0.2_173.T2D.glm
#
##for A in $(cat ${Output_dir}/prune/Obesity.1e-4.min0.2_173.T2D.txt)
##do
##grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid >> ${Output_dir}/prune/raw.Obesity.1e-4.min0.2_173.T2D.glm
##done
#
#head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
#cat ${Output_dir}/header ${Output_dir}/prune/raw.Obesity.1e-4.min0.2_173.T2D.glm > ${Output_dir}/prune/Obesity.1e-4.min0.2_173.T2D.glm
#rm -rf ${Output_dir}/prune/raw.Obesity.1e-4.min0.2_173.T2D.glm
#
#export CUTOFF="1e-4.min0.2_173.T2D"
#singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker_Obe_T2D.R
#
#echo ""
#echo "----------------------------"
#echo " Finish     1e-4.min0.2"
#echo "----------------------------"
#echo ""
