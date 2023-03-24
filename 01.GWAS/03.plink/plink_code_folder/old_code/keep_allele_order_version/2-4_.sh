#!/bin/sh
if [ $# -ne 1 ];then
        echo "Please enter Sample_Name"
               exit
fi
Sample=$1
GWAS_path="/scratch/x1997a11/GWAS/pdxen_AD/result_folder"
code_path="/scratch/x1997a11/GWAS/pdxen_AD/Code_folder"
Sample_folder="${GWAS_path}/QC_${Sample}"
#분석결과를 담을 파일
mkdir ${GWAS_path}/${Sample}_analysis_folder
analysis_folder="/${GWAS_path}/${Sample}_analysis_folder"
Case_pheno="${analysis_folder}/merge/Case_pheno.txt"
Control_pheno="${analysis_folder}/merge/Control_pheno.txt"
inversion="/scratch/x1997a11/GWAS/pdxen_AD/reference_folder/inversion.txt"
predixcan_folder="/scratch/x1997a11/GWAS/PrediXcan/MetaXcan-master/software"
Reference_folder="/scratch/x1997a11/REFERENCE"

#=================================================================


#snpsift dbsnp154
grep '^#' ${analysis_folder}/merge/logistic/anno/chr22_${Sample}_NoNA.vcf > ${analysis_folder}/merge/logistic/anno/merge_rmsnp_snpshift.vcf
for A in $(seq 1 22)
do
grep -v '^#' ${analysis_folder}/merge/logistic/anno/chr${A}_${Sample}_NoNA.vcf > ${analysis_folder}/merge/logistic/anno/nohead_chr${A}_${Sample}_NoNA.vcf
cat ${analysis_folder}/merge/logistic/anno/nohead_chr${A}_${Sample}_NoNA.vcf >> ${analysis_folder}/merge/logistic/anno/merge_rmsnp_snpshift.vcf
done




plink --bfile ${analysis_folder}/merge/${Sample}_NoNA \
 --freq case-control \
 --keep-allele-order \
 --out ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq

#anno_summary_file
#chr
awk '{print $1}'  ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq.frq.cc > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.chr
#posi
awk '{print $2}' ${analysis_folder}/merge/logistic/raw_${Sample}_NoNA_assoc.assoc.logistic > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.posi
#rsID
grep -v '^#' ${analysis_folder}/merge/logistic/anno/merge_rmsnp_snpshift.vcf | awk '{print $3}' > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.rsID
sed -i '1i\rsID' ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.rsID
#awk '{print $2}' ${analysis_folder}/merge/logistic/raw_${Sample}_NoNA_assoc.assoc.logistic > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.rsID
#A1-A2
awk '{print $3, $4}' ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq.frq.cc > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.A1-A2
#beta
awk '{print $7}' ${analysis_folder}/merge/logistic/raw_${Sample}_NoNA_assoc.assoc.logistic > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.beta
#pval
awk '{print $12}' ${analysis_folder}/merge/logistic/raw_${Sample}_NoNA_assoc.assoc.logistic > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.pval

echo '================================'
echo '                                '
echo '      Do making summary file    '
echo '                                '
echo '================================'


paste ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.rsID \
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq.frq.cc \
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.posi \
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.pval \
 > ${analysis_folder}/merge/logistic/anno/summary_result.csv
sed -i 's/\t/ /g' ${analysis_folder}/merge/logistic/anno/summary_result.csv
sed -i 's/\t/ /g' ${analysis_folder}/merge/logistic/anno/summary_result.csv
sed -i 's/        / /g' ${analysis_folder}/merge/logistic/anno/summary_result.csv
sed -i 's/      / /g' ${analysis_folder}/merge/logistic/anno/summary_result.csv
sed -i 's/    / /g' ${analysis_folder}/merge/logistic/anno/summary_result.csv
sed -i 's/  / /g' ${analysis_folder}/merge/logistic/anno/summary_result.csv


module purge
module load python/3.7
. /apps/applications/PYTHON/3.7/etc/profile.d/conda.sh
conda activate predixcan

#S-predixcan

mkdir ${analysis_folder}/merge/logistic/predixcan
echo '================================'
echo '                                '
echo '    Do making predixcan file    '
echo '                                '
echo '================================'
#Input file
paste ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.chr \
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.rsID \
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.A1-A2 \
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.beta \
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.pval \
 > ${analysis_folder}/merge/logistic/predixcan/${Sample}.predixcan

for A in $(seq 1 22)
do
grep -w "^${A}" ${analysis_folder}/merge/logistic/predixcan/${Sample}.predixcan > ${analysis_folder}/merge/logistic/predixcan/raw_chr${A}_${Sample}.predixcan
grep -v ";" ${analysis_folder}/merge/logistic/predixcan/raw_chr${A}_${Sample}.predixcan >  ${analysis_folder}/merge/logistic/predixcan/raw_rmsemi_chr${A}_${Sample}.predixcan
awk '{print $2, $3, $4, $5, $6}' ${analysis_folder}/merge/logistic/predixcan/raw_rmsemi_chr${A}_${Sample}.predixcan > ${analysis_folder}/merge/logistic/predixcan/raw_rmsemi_rmchr_chr${A}_${Sample}.predixcan
grep  '^rs' ${analysis_folder}/merge/logistic/predixcan/raw_rmsemi_rmchr_chr${A}_${Sample}.predixcan > ${analysis_folder}/merge/logistic/predixcan/chr${A}_${Sample}.predixcan
sed -i '1i\rsID\tA1\tA2\tBETA\tP' ${analysis_folder}/merge/logistic/predixcan/chr${A}_${Sample}.predixcan
gzip -f ${analysis_folder}/merge/logistic/predixcan/chr${A}_${Sample}.predixcan
rm ${analysis_folder}/merge/logistic/predixcan/raw*
done

python ${predixcan_folder}/SPrediXcan.py \
 --model_db_path ${Reference_folder}/prediXcan/weights/gtex_v7_Whole_Blood_imputed_europeans_tw_0.5_signif.db \
 --covariance ${Reference_folder}/prediXcan/covar_gtex/gtex_v7_Whole_Blood_imputed_eur_covariances.txt.gz \
 --gwas_folder ${analysis_folder}/merge/logistic/predixcan \
 --gwas_file_pattern ".*gz" \
 --snp_column rsID \
 --effect_allele_column A1 \
 --non_effect_allele_column A2 \
 --beta_column BETA \
 --pvalue_column P \
 --output_file ${analysis_folder}/merge/logistic/predixcan/result.csv
