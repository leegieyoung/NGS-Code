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
#=================================================================

#anno
mkdir ${analysis_folder}/merge/interaction

plink --bfile ${analysis_folder}/merge/${Sample}_NoNA \
 --extract ${GWAS_path}/19snp.list \
 --keep-allele-order \
 --make-bed \
 --out ${analysis_folder}/merge/interaction/19snp

#mds
plink --bfile ${analysis_folder}/merge/interaction/19snp \
 --logistic sex interaction \
 --covar ${analysis_folder}/merge/logistic.covar_mds_age_apoe.txt \
 --parameters 1-13, 25, 26 \
 --out ${analysis_folder}/merge/interaction/19snp_mds_sex_age_apoe

grep 'ADDxapoe' ${analysis_folder}/merge/interaction/19snp_mds_sex_age_apoe.assoc.logistic > ${analysis_folder}/merge/interaction/19snp_mds_sex_age_apoe.interaction

#dom
plink --bfile ${analysis_folder}/merge/interaction/19snp \
 --logistic sex dominant interaction \
 --covar ${analysis_folder}/merge/logistic.covar_mds_age_apoe.txt \
 --parameters 1-13, 25, 26 \
 --out ${analysis_folder}/merge/interaction/19snp_mds_sex_age_apoe_dom

grep 'DOMxapoe' ${analysis_folder}/merge/interaction/19snp_mds_sex_age_apoe_dom.assoc.logistic > ${analysis_folder}/merge/interaction/19snp_mds_sex_age_apoe_dom.interaction

#rec
plink --bfile ${analysis_folder}/merge/interaction/19snp \
 --logistic sex recessive interaction \
 --covar ${analysis_folder}/merge/logistic.covar_mds_age_apoe.txt \
 --parameters 1-13, 25, 26 \
 --out ${analysis_folder}/merge/interaction/19snp_mds_sex_age_apoe_rec

grep 'RECxapoe' ${analysis_folder}/merge/interaction/19snp_mds_sex_age_apoe_rec.assoc.logistic > ${analysis_folder}/merge/interaction/19snp_mds_sex_age_apoe.interaction

#Not mds
plink --bfile ${analysis_folder}/merge/interaction/19snp \
 --logistic sex interaction \
 --covar ${analysis_folder}/merge/logistic.covar_mds_age_apoe.txt \
 --parameters 1, 12-13, 25, 26 \
 --out ${analysis_folder}/merge/interaction/19snp_sex_age_apoe
grep 'ADDxapoe' ${analysis_folder}/merge/interaction/19snp_sex_age_apoe.assoc.logistic > ${analysis_folder}/merge/interaction/19snp_sex_age_apoe.interaction

#dom
plink --bfile ${analysis_folder}/merge/interaction/19snp \
 --logistic sex dominant interaction \
 --covar ${analysis_folder}/merge/logistic.covar_mds_age_apoe.txt \
 --parameters 1, 12-13, 25, 26 \
 --out ${analysis_folder}/merge/interaction/19snp_sex_age_apoe_dom

grep 'DOMxapoe' ${analysis_folder}/merge/interaction/19snp_sex_age_apoe_dom.assoc.logistic > ${analysis_folder}/merge/interaction/19snp_sex_age_apoe_dom.interaction

#rec
plink --bfile ${analysis_folder}/merge/interaction/19snp \
 --logistic sex recessive interaction \
 --covar ${analysis_folder}/merge/logistic.covar_mds_age_apoe.txt \
 --parameters 1, 12-13, 25, 26 \
 --out ${analysis_folder}/merge/interaction/19snp_sex_age_apoe_rec

grep 'RECxapoe' ${analysis_folder}/merge/interaction/19snp_sex_age_apoe_rec.assoc.logistic > ${analysis_folder}/merge/interaction/19snp_sex_age_apoe_rec.interaction

awk '$9 < 0.05 {print $0}' ${analysis_folder}/merge/interaction/*.interaction > ${analysis_folder}/merge/interaction/5e-2.pval

#logi
plink --bfile ${analysis_folder}/merge/interaction/19snp \
 --logistic sex \
 --covar ${analysis_folder}/merge/logistic.covar_mds_age_apoe.txt \
 --parameters 1-12, 14 \
--out ${analysis_folder}/merge/interaction/19snp_assoc

plink --bfile ${analysis_folder}/merge/interaction/19snp \
 --freq case-control \
 --keep-allele-order \
 --out ${analysis_folder}/merge/interaction/19snp_freq

mkdir ${analysis_folder}/merge/interaction/ppt
awk '{print $2}' ${analysis_folder}/merge/interaction/19snp_freq.frq.cc > ${analysis_folder}/merge/interaction/ppt/19snp.snp
awk '{print $4}' ${analysis_folder}/merge/interaction/19snp_freq.frq.cc > ${analysis_folder}/merge/interaction/ppt/19snp.A2
awk '{print $3}' ${analysis_folder}/merge/interaction/19snp_freq.frq.cc > ${analysis_folder}/merge/interaction/ppt/19snp.A1
awk '{print $7}' ${analysis_folder}/merge/interaction/19snp_freq.frq.cc > ${analysis_folder}/merge/interaction/ppt/19snp.Case
awk '{print $8}' ${analysis_folder}/merge/interaction/19snp_freq.frq.cc > ${analysis_folder}/merge/interaction/ppt/19snp.Control
awk '{print $5}' ${analysis_folder}/merge/interaction/19snp_freq.frq.cc > ${analysis_folder}/merge/interaction/ppt/19snp.Case_freq
awk '{print $6}' ${analysis_folder}/merge/interaction/19snp_freq.frq.cc > ${analysis_folder}/merge/interaction/ppt/19snp.Control_freq
grep 'ADD' ${analysis_folder}/merge/interaction/19snp_assoc.assoc.logistic > ${analysis_folder}/merge/interaction/19snp_assoc.assoc.ADD
sed -i '1i\chr\tsnp\tposi\tA1\tADD\twhat\twhat\twhat\tp' ${analysis_folder}/merge/interaction/19snp_assoc.assoc.ADD
awk '{print $9}' ${analysis_folder}/merge/interaction/19snp_assoc.assoc.ADD > ${analysis_folder}/merge/interaction/ppt/19snp.pval
awk '{print $9}' ${analysis_folder}/merge/interaction/19snp_mds_sex_age_apoe.interaction > ${analysis_folder}/merge/interaction/ppt/19snp_mds_sex_age_apoe.inter
sed -i '1i\mds_sex_age_apoe.pval' ${analysis_folder}/merge/interaction/ppt/19snp_mds_sex_age_apoe.inter
awk '{print $9}' ${analysis_folder}/merge/interaction/19snp_sex_age_apoe.interaction > ${analysis_folder}/merge/interaction/ppt/19snp_sex_age_apoe.inter
sed -i '1i\sex_age_apoe.pval' ${analysis_folder}/merge/interaction/ppt/19snp_sex_age_apoe.inter

paste -d '\t' ${analysis_folder}/merge/interaction/ppt/19snp.snp \
 ${analysis_folder}/merge/interaction/ppt/19snp.Case \
 ${analysis_folder}/merge/interaction/ppt/19snp.Control \
 ${analysis_folder}/merge/interaction/ppt/19snp.A2 \
 ${analysis_folder}/merge/interaction/ppt/19snp.A1 \
 ${analysis_folder}/merge/interaction/ppt/19snp.Case_freq \
 ${analysis_folder}/merge/interaction/ppt/19snp.Control_freq \
 ${analysis_folder}/merge/interaction/ppt/19snp.pval \
 ${analysis_folder}/merge/interaction/ppt/19snp_mds_sex_age_apoe.inter \
 ${analysis_folder}/merge/interaction/ppt/19snp_sex_age_apoe.inter > ${analysis_folder}/merge/interaction/ppt/result.csv




