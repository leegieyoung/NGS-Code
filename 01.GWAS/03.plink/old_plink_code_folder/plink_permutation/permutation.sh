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
echo "${code_path}/2_Single_version" > ${analysis_folder}/merge/2_single_version

plink --bfile ${QC_folder}/MaMi/QC_${Sample} \
 --keep-allele-order \
 --assoc \
 --out ${analysis_folder}/merge/raw_${Sample}_assoc

#grep -v "NA" ${analysis_folder}/merge/raw_${Sample}_assoc.assoc | grep -v "CHR" | awk '{print $2}' > ${analysis_folder}/merge/extract_SNP_list.txt
awk '$9 != "NA" {print $0}' ${analysis_folder}/merge/raw_${Sample}_assoc.assoc | grep -v "CHR" | awk '{print $2}' > ${analysis_folder}/merge/extract_SNP_list.txt


#==========================================================

plink --bfile ${QC_folder}/MaMi/QC_${Sample} \
 --extract ${analysis_folder}/merge/extract_SNP_list.txt \
 --keep-allele-order \
 --make-bed \
 --out ${analysis_folder}/merge/raw_${Sample}_NoNA
#==========================================================


#prune
plink --bfile ${analysis_folder}/merge/raw_${Sample}_NoNA \
 --keep-allele-order \
 --exclude ${inversion} \
 --range \
 --indep-pairwise 50 5 0.2 \
 --out ${analysis_folder}/merge/raw_indepSNP

#assoc

#genome
plink --bfile ${analysis_folder}/merge/raw_${Sample}_NoNA \
 --extract ${analysis_folder}/merge/raw_indepSNP.prune.in \
 --keep-allele-order \
 --genome \
 --out ${analysis_folder}/merge/raw_${Sample}_NoNA_genome

#MDS
mkdir ${analysis_folder}/merge/MDS
plink --bfile ${analysis_folder}/merge/raw_${Sample}_NoNA \
 --extract ${analysis_folder}/merge/raw_indepSNP.prune.in \
 --keep-allele-order \
 --make-bed \
 --out ${analysis_folder}/merge/raw_indep_${Sample}_NoNA

awk '{print $1, $2, $6}' ${analysis_folder}/merge/raw_${Sample}_NoNA.fam > ${analysis_folder}/merge/raw_pheno.txt
awk '$3 > 1 {print $0}' ${analysis_folder}/merge/raw_pheno.txt > ${analysis_folder}/merge/Case_pheno.txt
awk '$3 < 2 && $3 > 0 {print $0}' ${analysis_folder}/merge/raw_pheno.txt > ${analysis_folder}/merge/Control_pheno.txt 

plink --bfile ${analysis_folder}/merge/raw_indep_${Sample}_NoNA \
 --read-genome ${analysis_folder}/merge/raw_${Sample}_NoNA_genome.genome \
 --cluster --mds-plot 10 \
 --keep-allele-order \
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
 --set-missing-var-ids @:# \
 --keep-allele-order \
 --out ${analysis_folder}/merge/PCA/${Sample}_NoNA_PCA

awk '{print $1, $3, $4, $5, $6}' ${analysis_folder}/merge/PCA/${Sample}_NoNA_PCA.eigenvec > ${analysis_folder}/merge/PCA/raw_${Sample}_PCA.csv
sed -i '1i\name PC1 PC2 PC3 PC4' ${analysis_folder}/merge/PCA/raw_${Sample}_PCA.csv
paste -d '\t' ${analysis_folder}/merge/PCA/raw_${Sample}_PCA.csv ${analysis_folder}/merge/PCA/${Sample}_NoNA_PCA.eigenval > ${analysis_folder}/merge/PCA/raw_${Sample}_PCA-eigenval.csv

#Manhattan - QQplot code
cp ${Rcode}/Manhattan_plot.R ${analysis_folder}/merge/logistic/Manhattan_plot.R


#logistic regression (Odd Ratio가 아닌 Beta)
mkdir ${analysis_folder}/merge/logistic
awk '{print $1, $2, $4, $5, $6, $7 ,$8 ,$9 ,$10 ,$11, $12, $13}' ${analysis_folder}/merge/MDS/${Sample}_NoNA_genome_MDS.mds > ${analysis_folder}/merge/logistic.covar_mds.txt

plink --bfile ${analysis_folder}/merge/raw_${Sample}_NoNA \
 --covar ${analysis_folder}/merge/logistic.covar_mds.txt \
 --logistic beta\
 --sex \
 --keep-allele-order \
 --hide-covar \
 --ci 0.95 \
 --out ${analysis_folder}/merge/logistic/raw_${Sample}_NoNA_assoc

awk '!/'NA'/' ${analysis_folder}/merge/logistic/raw_${Sample}_NoNA_assoc.assoc.logistic > ${analysis_folder}/merge/logistic/${Sample}_NoNA_assoc.assoc.logistic

#haploview folder
mkdir ${analysis_folder}/merge/haploview

sed -i 's/        / /g' ${analysis_folder}/merge/logistic/${Sample}_NoNA_assoc.assoc.logistic
sed -i 's/    / /g' ${analysis_folder}/merge/logistic/${Sample}_NoNA_assoc.assoc.logistic
sed -i 's/  / /g' ${analysis_folder}/merge/logistic/${Sample}_NoNA_assoc.assoc.logistic


awk '$12 < 0.05 {print $0}' ${analysis_folder}/merge/logistic/${Sample}_NoNA_assoc.assoc.logistic > ${analysis_folder}/merge/logistic/low_${Sample}_NoNA_assoc.assoc.logistic
sed -i '1i\ CHR SNP BP A1 TEST NMISS BETA SE L95 U95 STAT P' ${analysis_folder}/merge/logistic/low_${Sample}_NoNA_assoc.assoc.logistic

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
 --keep-allele-order \
 --out ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq

plink --bfile ${analysis_folder}/merge/${Sample}_NoNA \
 --freq \
 --keep-allele-order \
 --out ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq
 
sed -i 's/                                                                     / /g' ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq.frq.cc
sed -i 's/      / /g' ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq.frq.cc
sed -i 's/        / /g' ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq.frq.cc
sed -i 's/\t/ /g' ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq.frq.cc
sed -i 's/   / /g' ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq.frq.cc
sed -i 's/  / /g' ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq.frq.cc

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
#
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

gzip ${analysis_folder}/merge/logistic/anno/${Sample}.txt
