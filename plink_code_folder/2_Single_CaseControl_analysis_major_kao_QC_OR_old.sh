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
 --assoc \
 --threads 64 \
 --out ${analysis_folder}/merge/raw_${Sample}_assoc

#grep -v "NA" ${analysis_folder}/merge/raw_${Sample}_assoc.assoc | grep -v "CHR" | awk '{print $2}' > ${analysis_folder}/merge/extract_SNP_list.txt
awk '$9 != "NA" {print $0}' ${analysis_folder}/merge/raw_${Sample}_assoc.assoc | grep -v "CHR" | awk '{print $2}' > ${analysis_folder}/merge/extract_SNP_list.txt


#==========================================================

plink --bfile ${QC_folder}/MaMi/QC_${Sample} \
 --extract ${analysis_folder}/merge/extract_SNP_list.txt \
 --make-bed \
 --threads 64 \
 --out ${analysis_folder}/merge/raw_snpQC_${Sample}_NoNA
#==========================================================


#prune
plink --bfile ${analysis_folder}/merge/raw_snpQC_${Sample}_NoNA \
 --exclude ${inversion} \
 --range \
 --indep-pairwise 50 5 0.2 \
 --threads 64 \
 --out ${analysis_folder}/merge/raw_snpQC_indepSNP

#IBD
plink --bfile ${analysis_folder}/merge/raw_snpQC_${Sample}_NoNA \
 --extract ${analysis_folder}/merge/raw_snpQC_indepSNP.prune.in \
 --genome \
 --min 0.2 \
 --threads 64 \
 --out ${analysis_folder}/merge/raw_remove

#Remove list
awk '{print $1}' ${analysis_folder}/merge/raw_remove.genome | sed '1d' > ${analysis_folder}/merge/raw_remove.txt
paste ${analysis_folder}/merge/raw_remove.txt ${analysis_folder}/merge/raw_remove.txt > ${analysis_folder}/merge/raw_remove.list

#Remove sample
plink --bfile ${analysis_folder}/merge/raw_snpQC_${Sample}_NoNA \
 --remove ${analysis_folder}/merge/raw_remove.list \
 --make-bed \
 --threads 64 \
 --out ${analysis_folder}/merge/raw_snpQC_IBDQC_${Sample}_NoNA

#Remove Variant in R2 single
plink --bfile ${analysis_folder}/merge/raw_snpQC_IBDQC_${Sample}_NoNA \
  --show-tags ${analysis_folder}/merge/raw_snpQC_indepSNP.prune.in \
  --list-all \
 --tag-kb 250 \
 --tag-r2 0.8 \
 --threads 64 \
 --out ${analysis_folder}/merge/raw_tagQC

awk '{print $1, $8}' ${analysis_folder}/merge/raw_tagQC.tags.list | grep 'NONE' | awk '{print $1}' > ${analysis_folder}/merge/raw_exclude.list

plink --bfile ${analysis_folder}/merge/raw_snpQC_IBDQC_${Sample}_NoNA \
 --exclude ${analysis_folder}/merge/raw_exclude.list \
 --make-bed \
 --threads 64 \
 --out ${analysis_folder}/merge/raw_${Sample}_NoNA

#prune for MDS and PCA
plink --bfile ${analysis_folder}/merge/raw_${Sample}_NoNA \
 --threads 64 \
 --exclude ${inversion} \
 --range \
 --indep-pairwise 50 5 0.2 \
 --out ${analysis_folder}/merge/raw_indepSNP

plink --bfile ${analysis_folder}/merge/raw_${Sample}_NoNA \
 --extract ${analysis_folder}/merge/raw_indepSNP.prune.in \
 --threads 64 \
 --genome \
 --out ${analysis_folder}/merge/raw_${Sample}_NoNA_genome

#MDS
mkdir ${analysis_folder}/merge/MDS
plink --bfile ${analysis_folder}/merge/raw_${Sample}_NoNA \
 --extract ${analysis_folder}/merge/raw_indepSNP.prune.in \
 --threads 64 \
 --make-bed \
 --out ${analysis_folder}/merge/raw_indep_${Sample}_NoNA

awk '{print $1, $2, $6}' ${analysis_folder}/merge/raw_${Sample}_NoNA.fam > ${analysis_folder}/merge/raw_pheno.txt
awk '$3 > 1 {print $0}' ${analysis_folder}/merge/raw_pheno.txt > ${analysis_folder}/merge/Case_pheno.txt
awk '$3 < 2 && $3 > 0 {print $0}' ${analysis_folder}/merge/raw_pheno.txt > ${analysis_folder}/merge/Control_pheno.txt 

plink --bfile ${analysis_folder}/merge/raw_indep_${Sample}_NoNA \
 --read-genome ${analysis_folder}/merge/raw_${Sample}_NoNA_genome.genome \
 --cluster --mds-plot 10 \
 --threads 64 \
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
 --threads 64 \
 --out ${analysis_folder}/merge/PCA/${Sample}_NoNA_PCA

awk '{print $1, $3, $4, $5, $6}' ${analysis_folder}/merge/PCA/${Sample}_NoNA_PCA.eigenvec > ${analysis_folder}/merge/PCA/raw_${Sample}_PCA.csv
sed -i '1i\name PC1 PC2 PC3 PC4' ${analysis_folder}/merge/PCA/raw_${Sample}_PCA.csv
paste -d '\t' ${analysis_folder}/merge/PCA/raw_${Sample}_PCA.csv ${analysis_folder}/merge/PCA/${Sample}_NoNA_PCA.eigenval > ${analysis_folder}/merge/PCA/raw_${Sample}_PCA-eigenval.csv

awk '{print $6}' ${analysis_folder}/merge/raw_${Sample}_NoNA.fam > ${analysis_folder}/merge/PCA/${Sample}_NoNA.pheno
sed -i -e 's/2/Case/g' -e 's/1/Control/g' ${analysis_folder}/merge/PCA/${Sample}_NoNA.pheno
sed -i '1iname' ${analysis_folder}/merge/PCA/${Sample}_NoNA.pheno
awk '{print $1}' ${analysis_folder}/merge/raw_${Sample}_NoNA.fam > ${analysis_folder}/merge/PCA/${Sample}_NoNA.sample
sed -i '1isample' ${analysis_folder}/merge/PCA/${Sample}_NoNA.sample
awk '{print $3, $4}' ${analysis_folder}/merge/PCA/${Sample}_NoNA_PCA.eigenvec > ${analysis_folder}/merge/PCA/raw_${Sample}_PCA.PC12
sed -i '1iPC1 PC2' ${analysis_folder}/merge/PCA/raw_${Sample}_PCA.PC12
paste -d ' ' ${analysis_folder}/merge/PCA/${Sample}_NoNA.pheno ${analysis_folder}/merge/PCA/raw_${Sample}_PCA.PC12 ${analysis_folder}/merge/PCA/${Sample}_NoNA.sample > ${analysis_folder}/merge/PCA/PCA.txt

cp ${Rcode}/pca.R ${analysis_folder}/merge/PCA/


#logistic regression (Odd Ratio가 아닌 Beta)
mkdir ${analysis_folder}/merge/logistic
awk '{print $1, $2, $4, $5, $6, $7 ,$8 ,$9 ,$10 ,$11, $12, $13}' ${analysis_folder}/merge/MDS/${Sample}_NoNA_genome_MDS.mds > ${analysis_folder}/merge/logistic.covar_mds.txt

plink --bfile ${analysis_folder}/merge/raw_${Sample}_NoNA \
 --covar ${analysis_folder}/merge/logistic.covar_mds.txt \
 --logistic \
 --sex \
 --hide-covar \
 --ci 0.95 \
 --threads 64 \
 --out ${analysis_folder}/merge/logistic/raw_${Sample}_NoNA_assoc

#awk '!/'NA'/' ${analysis_folder}/merge/logistic/raw_${Sample}_NoNA_assoc.assoc.logistic > ${analysis_folder}/merge/logistic/${Sample}_NoNA_assoc.assoc.logistic
awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' ${analysis_folder}/merge/logistic/raw_${Sample}_NoNA_assoc.assoc.logistic > ${analysis_folder}/merge/logistic/${Sample}_NoNA_assoc.assoc.logistic

#Manhattan - QQplot code
cp ${Rcode}/Manhattan_plot.R ${analysis_folder}/merge/logistic/

#haploview folder
mkdir ${analysis_folder}/merge/haploview

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
 --threads 64 \
 --out ${analysis_folder}/merge/logistic/anno/raw_${Sample}_NoNA_freq

plink --bfile ${analysis_folder}/merge/${Sample}_NoNA \
 --freq \
 --threads 64 \
 --out ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq
 
awk '{print $1, $2, $3, $4, $5, $6, $7, $8}' ${analysis_folder}/merge/logistic/anno/raw_${Sample}_NoNA_freq.frq.cc > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq.frq.cc

#anno_summary_file
#chr
awk '{print $1}'  ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq.frq.cc > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.chr
#posi
awk '{print $3}' ${analysis_folder}/merge/logistic/raw_${Sample}_NoNA_assoc.assoc.logistic > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.posi
#rsID
awk '{print $2}' ${analysis_folder}/merge/logistic/raw_${Sample}_NoNA_assoc.assoc.logistic > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.rsID
#A1-A2
awk '{print $3, $4}' ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq.frq.cc > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.A1-A2
awk '{print $3}' ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq.frq.cc > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.A1
awk '{print $4}' ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA_freq.frq.cc > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.A2
#OR
awk '{print $7}' ${analysis_folder}/merge/logistic/raw_${Sample}_NoNA_assoc.assoc.logistic > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.OR
#pval
awk '{print $12}' ${analysis_folder}/merge/logistic/raw_${Sample}_NoNA_assoc.assoc.logistic > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.pval
#SE
awk '{print $8}' ${analysis_folder}/merge/logistic/raw_${Sample}_NoNA_assoc.assoc.logistic > ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.SE
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
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.OR \
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
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.OR \
 ${analysis_folder}/merge/logistic/anno/${Sample}_NoNA.SE > ${analysis_folder}/merge/logistic/FUMA/fuma_input.txt

sed -i '1d' ${analysis_folder}/merge/logistic/FUMA/fuma_input.txt
sed -i '1ichromosome\tposition\tSNP\tP-value\tA1\tA2\tOR\tSE' ${analysis_folder}/merge/logistic/FUMA/fuma_input.txt

echo ""
echo "===========================Delete Used Data==========================="
echo ""
rm ${analysis_folder}/merge/logistic/anno/*_freq.frq
rm ${analysis_folder}/merge/logistic/anno/*_freq.frq.cc
rm ${analysis_folder}/merge/logistic/anno/${Sample}.txt
rm ${analysis_folder}/merge/logistic/raw_${Sample}_NoNA_assoc.assoc.logistic
rm ${analysis_folder}/merge/raw_file/raw_${Sample}_assoc.assoc

