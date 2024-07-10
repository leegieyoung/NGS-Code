#!/bin/sh
if [ $# -ne 1 ];then
        echo "Please enter Sample_Name"
               exit
fi
Sample=$1
GWAS_path="/scratch/hpc46a05/GWAS/result"
code_path="/scratch/hpc46a05/GWAS/Code/plink_code_folder"
QC_folder="${GWAS_path}/QC_${Sample}"
#imputed_sample folder
beagle_folder="/scratch/hpc46a05/GWAS/beagle_result"
Sample_folder="${beagle_folder}/3.${Sample}_imputed_${Sample}"
#분석결과를 담는 파일
mkdir ${GWAS_path}/${Sample}_analysis_folder
analysis_folder="/${GWAS_path}/${Sample}_analysis_folder"
Case_pheno="${analysis_folder}/merge/Case_pheno.txt"
Control_pheno="${analysis_folder}/merge/Control_pheno.txt"
inversion="/scratch/hpc46a05/REFERENCE/inversion.txt"
Rcode="/scratch/hpc46a05/GWAS/Rcode"
#=================================================================

mkdir ${analysis_folder}/merge
#=========What is Code ? =======================================
#anno

echo '================================'
echo '                                '
echo '    Do making FUMA input file '
echo '                                '
echo '================================'

mkdir ${analysis_folder}/merge/logistic/FUMA
paste -d '\t' ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.chr \
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.posi \
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.rsID \
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.pval \
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.A1 \
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.A2 \
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.beta \
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.SE \
 > ${analysis_folder}/merge/logistic/FUMA/fuma_input.txt

sed -i '1d' ${analysis_folder}/merge/logistic/FUMA/fuma_input.txt
sed -i '1ichromosome\tposition\tSNP\tP-value\tA1\tA2\tBETA\tSE' ${analysis_folder}/merge/logistic/FUMA/fuma_input.txt

echo ""
echo "===========================Delete Used Data==========================="
echo ""
rm ${analysis_folder}/merge/logistic/anno/*_freq.frq
rm ${analysis_folder}/merge/logistic/anno/*_freq.frq.cc
rm ${analysis_folder}/merge/logistic/anno/${Sample}.txt
rm ${analysis_folder}/merge/logistic/${Sample}_NoNA_assoc.assoc.logistic
rm ${analysis_folder}/merge/raw_file/raw_${Sample}_assoc.assoc

 
