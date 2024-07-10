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
mkdir ${analysis_folder}/merge/logistic/anno

plink --bfile ${analysis_folder}/merge/${Sample}_NoNA \
 --freq case-control \
 --keep-allele-order \
 --out ${analysis_folder}/merge/logistic/anno/raw_${Sample}_NoNA_freq

plink --bfile ${analysis_folder}/merge/${Sample}_NoNA \
 --freq \
 --keep-allele-order \
 --out ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq

awk '{print $1, $2, $3, $4, $5, $6, $7, $8}' ${analysis_folder}/merge/logistic/anno/raw_${Sample}_NoNA_freq.frq.cc > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq.frq.cc

#anno_summary_file
#chr
awk '{print $1}'  ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq.frq.cc > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.chr
#posi
awk '{print $3}' ${analysis_folder}/merge/logistic/test_mds_age_sex.assoc.logistic > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.posi
#rsID
awk '{print $2}' ${analysis_folder}/merge/logistic/test_mds_age_sex.assoc.logistic > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.rsID
#A1-A2
awk '{print $3, $4}' ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq.frq.cc > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.A1-A2
awk '{print $3}' ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq.frq.cc > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.A1
awk '{print $4}' ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq.frq.cc > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.A2
#beta
awk '{print $7}' ${analysis_folder}/merge/logistic/test_mds_age_sex.assoc.logistic > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.beta
#pval
awk '{print $12}' ${analysis_folder}/merge/logistic/test_mds_age_sex.assoc.logistic > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.pval
#SE
awk '{print $8}' ${analysis_folder}/merge/logistic/test_mds_age_sex.assoc.logistic > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.SE
#maf
awk '{print $5}' ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq.frq > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.maf


echo '================================'
echo '                                '
echo '      Do making summary file    '
echo '                                '
echo '================================'

paste ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq.frq.cc \
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.posi \
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.pval \
 > ${analysis_folder}/merge/logistic/anno/summary_result.csv

echo '================================'
echo '                                '
echo ' Do making predixcan input file '
echo '                                '
echo '================================'
paste -d '\t' ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.rsID \
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.chr \
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.posi \
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.A1-A2 \
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.maf \
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.beta \
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.pval \
 > ${analysis_folder}/merge/logistic/anno/raw_${Sample}.txt

grep -v 'NA' ${analysis_folder}/merge/logistic/anno/raw_${Sample}.txt > ${analysis_folder}/merge/logistic/anno/${Sample}.txt
sed -i 's/ /\t/g' ${analysis_folder}/merge/logistic/anno/${Sample}.txt

yes n | gzip ${analysis_folder}/merge/logistic/anno/${Sample}.txt 

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
rm ${analysis_folder}/merge/raw_file/raw_${Sample}_assoc.assoc

 
