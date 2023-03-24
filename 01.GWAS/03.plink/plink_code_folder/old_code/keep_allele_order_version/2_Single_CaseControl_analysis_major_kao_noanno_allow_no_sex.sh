#!/bin/sh
if [ $# -ne 1 ];then
        echo "Please enter Sample_Name"
               exit
fi
Sample=$1
GWAS_path="/data/keeyoung/GWAS/result"
code_path="/data/keeyoung/GWAS/Code/plink_code_folder/keep_allele_order_version"
QC_folder="${GWAS_path}/QC_${Sample}"
#imputed_sample folder
beagle_folder="/data/keeyoung/GWAS/beagle_result"
Sample_folder="${beagle_folder}/3.${Sample}_imputed_${Sample}"
#분석결과를 담는 파일
mkdir ${GWAS_path}/${Sample}_analysis_folder
analysis_folder="/${GWAS_path}/${Sample}_analysis_folder"
Case_pheno="${analysis_folder}/merge/Case_pheno.txt"
Control_pheno="${analysis_folder}/merge/Control_pheno.txt"
inversion="/data/keeyoung/REFERENCE/inversion.txt"
#=================================================================

mkdir ${analysis_folder}/merge
#=========What is Code ? =======================================
echo "${code_path}/3_Single_CaseControl_analysis_major_kao.sh" > ${analysis_folder}/merge/Used_2_Single_CaseControl_analysis_major_kao.sh

plink --bfile ${QC_folder}/MaMi/QC_${Sample} \
 --keep-allele-order \
 --assoc \
 --allow-no-sex \
 --out ${analysis_folder}/merge/raw_${Sample}_assoc

#grep -v "NA" ${analysis_folder}/merge/raw_${Sample}_assoc.assoc | grep -v "CHR" | awk '{print $2}' > ${analysis_folder}/merge/extract_SNP_list.txt
awk '$9 != "NA" {print $0}' ${analysis_folder}/merge/raw_${Sample}_assoc.assoc | grep -v "CHR" | awk '{print $2}' > ${analysis_folder}/merge/extract_SNP_list.txt


#==========================================================

plink --bfile ${QC_folder}/MaMi/QC_${Sample} \
 --extract ${analysis_folder}/merge/extract_SNP_list.txt \
 --keep-allele-order \
 --make-bed \
 --maf 0.01 \
 --allow-no-sex \
 --out ${analysis_folder}/merge/raw_${Sample}_NoNA
#==========================================================


#prune
plink --bfile ${analysis_folder}/merge/raw_${Sample}_NoNA \
 --keep-allele-order \
 --exclude ${inversion} \
 --range \
 --allow-no-sex \
 --indep-pairwise 50 5 0.2 \
 --out ${analysis_folder}/merge/raw_indepSNP

#assoc
plink --bfile ${analysis_folder}/merge/raw_${Sample}_NoNA \
 --keep-allele-order \
 --assoc \
 --allow-no-sex \
 --out ${analysis_folder}/merge/raw_${Sample}_NoNA_assoc

awk '!/'NA'/' ${analysis_folder}/merge/raw_${Sample}_NoNA_assoc.assoc > ${analysis_folder}/merge/${Sample}_NoNA_assoc.assoc

#genome
plink --bfile ${analysis_folder}/merge/raw_${Sample}_NoNA \
 --extract ${analysis_folder}/merge/raw_indepSNP.prune.in \
 --keep-allele-order \
 --genome \
 --allow-no-sex \
 --out ${analysis_folder}/merge/raw_${Sample}_NoNA_genome

#MDS
mkdir ${analysis_folder}/merge/MDS
plink --bfile ${analysis_folder}/merge/raw_${Sample}_NoNA \
 --extract ${analysis_folder}/merge/raw_indepSNP.prune.in \
 --keep-allele-order \
 --make-bed \
 --allow-no-sex \
 --out ${analysis_folder}/merge/raw_indep_${Sample}_NoNA

awk '{print $1, $2, $6}' ${analysis_folder}/merge/raw_${Sample}_NoNA.fam > ${analysis_folder}/merge/raw_pheno.txt
awk '$3 > 1 {print $0}' ${analysis_folder}/merge/raw_pheno.txt > ${analysis_folder}/merge/Case_pheno.txt
awk '$3 < 2 && $3 > 0 {print $0}' ${analysis_folder}/merge/raw_pheno.txt > ${analysis_folder}/merge/Control_pheno.txt 

plink --bfile ${analysis_folder}/merge/raw_indep_${Sample}_NoNA \
 --read-genome ${analysis_folder}/merge/raw_${Sample}_NoNA_genome.genome \
 --cluster --mds-plot 10 \
 --keep-allele-order \
 --allow-no-sex \
 --out ${analysis_folder}/merge/MDS/${Sample}_NoNA_genome_MDS

cat ${Case_pheno} ${Control_pheno} > ${analysis_folder}/merge/MDS/${Sample}_pheno.txt
awk '{print $1, $2, "Control"}' ${analysis_folder}/merge/MDS/${Sample}_pheno.txt > ${analysis_folder}/merge/MDS/Control_pheno.txt
awk '{print $1, $2, "Case"}' ${analysis_folder}/merge/MDS/${Sample}_pheno.txt > ${analysis_folder}/merge/MDS/Case_pheno.txt
cat ${analysis_folder}/merge/MDS/Control_pheno.txt ${analysis_folder}/merge/MDS/Case_pheno.txt | sed -e '1i\FID IID pheno' > ${analysis_folder}/merge/MDS/phenofile.txt

#PCA
mkdir ${analysis_folder}/merge/PCA
plink --bfile ${analysis_folder}/merge/raw_${Sample}_NoNA \
 --double-id \
 --pca 10 \
 --allow-no-sex \
 --set-missing-var-ids @:# \
 --keep-allele-order \
 --out ${analysis_folder}/merge/PCA/${Sample}_NoNA_PCA

awk '{print $1, $3, $4, $5, $6}' ${analysis_folder}/merge/PCA/${Sample}_NoNA_PCA.eigenvec > ${analysis_folder}/merge/PCA/raw_${Sample}_PCA.csv
sed -i '1i\name PC1 PC2 PC3 PC4' ${analysis_folder}/merge/PCA/raw_${Sample}_PCA.csv

#Manhattan - QQplot code
cp /data/keeyoung/Rcode/Manhattan_plot.R ${analysis_folder}/merge/

#logistic regression (Odd Ratio가 아닌 Beta)
mkdir ${analysis_folder}/merge/logistic
awk '{print $1, $2, $4, $5, $6, $7 ,$8 ,$9 ,$10 ,$11, $12, $13}' ${analysis_folder}/merge/MDS/${Sample}_NoNA_genome_MDS.mds > ${analysis_folder}/merge/logistic.convar_mds.txt

plink --bfile ${analysis_folder}/merge/raw_${Sample}_NoNA \
 --covar ${analysis_folder}/merge/logistic.convar_mds.txt \
 --logistic beta\
 --allow-no-sex \
 --keep-allele-order \
 --hide-covar \
 --ci 0.95 \
 --out ${analysis_folder}/merge/logistic/raw_${Sample}_NoNA_assoc

awk '!/'NA'/' ${analysis_folder}/merge/logistic/raw_${Sample}_NoNA_assoc.assoc.logistic > ${analysis_folder}/merge/logistic/${Sample}_NoNA_assoc.assoc.logistic

#haploview folder
mkdir ${analysis_folder}/merge/haploview

sed -i 's/        / /g' ${analysis_folder}/merge/logistic/${Sample}_NoNA_assoc.assoc.logistic
sed -i 's/    / /g' ${analysis_folder}/merge/${Sample}_NoNA_assoc.assoc
sed -i 's/    / /g' ${analysis_folder}/merge/logistic/${Sample}_NoNA_assoc.assoc.logistic
sed -i 's/  / /g' ${analysis_folder}/merge/logistic/${Sample}_NoNA_assoc.assoc.logistic
sed -i 's/  / /g' ${analysis_folder}/merge/${Sample}_NoNA_assoc.assoc


awk '$12 < 0.05 {print $0}' ${analysis_folder}/merge/logistic/${Sample}_NoNA_assoc.assoc.logistic > ${analysis_folder}/merge/logistic/low_${Sample}_NoNA_assoc.assoc.logistic
sed -i '1i\ CHR SNP BP A1 TEST NMISS BETA SE L95 U95 STAT P' ${analysis_folder}/merge/logistic/low_${Sample}_NoNA_assoc.assoc.logistic
awk '$9 < 0.05 {print $0}' ${analysis_folder}/merge/${Sample}_NoNA_assoc.assoc > ${analysis_folder}/merge/low_${Sample}_NoNA_assoc.assoc
sed -i '1i\ CHR SNP BP A1 F_A F_U A2 CHISQ P OR' ${analysis_folder}/merge/low_${Sample}_NoNA_assoc.assoc


mv ${analysis_folder}/merge/raw_${Sample}_NoNA.bed ${analysis_folder}/merge/${Sample}_NoNA.bed
mv ${analysis_folder}/merge/raw_${Sample}_NoNA.bim ${analysis_folder}/merge/${Sample}_NoNA.bim
mv ${analysis_folder}/merge/raw_${Sample}_NoNA.fam ${analysis_folder}/merge/${Sample}_NoNA.fam
mv ${analysis_folder}/merge/raw_${Sample}_NoNA.log ${analysis_folder}/merge/${Sample}_NoNA.log
#rm ${analysis_folder}/merge/raw*
mkdir ${analysis_folder}/merge/raw_file
mv ${analysis_folder}/merge/raw* ${analysis_folder}/merge/raw_file

#anno
mkdir ${analysis_folder}/merge/logistic/anno

plink --bfile ${analysis_folder}/merge/${Sample}_NoNA \
 --freq case-control \
 --allow-no-sex \
 --keep-allele-order \
 --out ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq

#anno_summary_file
#chr
awk '{print $1}'  ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq.frq.cc > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.chr
#posi
awk '{print $3}' ${analysis_folder}/merge/logistic/raw_${Sample}_NoNA_assoc.assoc.logistic > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.posi
#rsID
awk '{print $2}' ${analysis_folder}/merge/logistic/raw_${Sample}_NoNA_assoc.assoc.logistic > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.rsID
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

paste ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq.frq.cc \
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.posi \
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.pval \
 > ${analysis_folder}/merge/logistic/anno/summary_result.csv


