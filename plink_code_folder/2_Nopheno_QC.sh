#!/bin/sh
if [ $# -ne 2 ];then
        echo "Please enter Sample_Name, thread(s)"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
#beagle 진행시Same position 제거를 안했을 경우 사용하는 코드
Sample=$1
Thread=$2
Dir="/mnt/nas/gylee/0.GWAS"
impute_dir="${Dir}/1.Input/${Sample}/0.Impute"
QC_dir="${Dir}/1.Input/${Sample}/1.QC"
MaMi_dir="${Dir}/1.Input/${Sample}/2.MaMi"
inversion="${Dir}/REFERENCE/inversion.txt"
Rcode="/mnt/nas/gylee/0.GWAS/Code/R"
#=================================================================

mkdir -p ${Dir}/2.plink_result/${Sample}/temp
mkdir -p ${Dir}/2.plink_result/${Sample}/result
ana_dir="${Dir}/2.plink_result/${Sample}"


#=========What is Code ? =======================================
echo "${code_path}/2_Single_version" > ${ana_dir}/2_single_version

plink --bfile ${MaMi_dir}/QC_Imputed_${Sample} \
 --assoc \
 --threads ${Thread} \
 --keep-allele-order \
 --allow-no-sex \
 --out ${ana_dir}/temp/1.raw_${Sample}_assoc

#grep -v "NA" ${ana_dir}/raw_${Sample}_assoc.assoc | grep -v "CHR" | awk '{print $2}' > ${ana_dir}/extract_SNP_list.txt
awk '$9 != "NA" {print $0}' ${ana_dir}/temp/1.raw_${Sample}_assoc.assoc | grep -v "CHR" | awk '{print $2}' > ${ana_dir}/extract_SNP_list.txt


#==========================================================

plink --bfile ${MaMi_dir}/QC_Imputed_${Sample} \
 --extract ${ana_dir}/extract_SNP_list.txt \
 --make-bed \
 --threads ${Thread} \
 --keep-allele-order \
 --allow-no-sex \
 --out ${ana_dir}/temp/2.raw_snpQC_${Sample}_NoNA
#==========================================================


#prune
plink --bfile ${ana_dir}/temp/2.raw_snpQC_${Sample}_NoNA \
 --exclude ${inversion} \
 --range \
 --indep-pairwise 50 5 0.2 \
 --threads ${Thread} \
 --out ${ana_dir}/temp/3.raw_snpQC_indepSNP

#IBD
plink --bfile ${ana_dir}/temp/2.raw_snpQC_${Sample}_NoNA \
 --extract ${ana_dir}/temp/3.raw_snpQC_indepSNP.prune.in \
 --genome \
 --min 0.2 \
 --threads ${Thread} \
 --allow-no-sex \
 --out ${ana_dir}/temp/4.raw_remove

#Remove list
awk '{print $1}' ${ana_dir}/temp/4.raw_remove.genome | sed '1d' > ${ana_dir}/temp/4.raw_remove.txt
paste ${ana_dir}/temp/4.raw_remove.txt ${ana_dir}/temp/4.raw_remove.txt > ${ana_dir}/temp/4.raw_remove.list

#Remove sample
plink --bfile ${ana_dir}/temp/2.raw_snpQC_${Sample}_NoNA \
 --remove ${ana_dir}/temp/4.raw_remove.list \
 --make-bed \
 --threads ${Thread} \
 --allow-no-sex \
 --out ${ana_dir}/temp/5.raw_snpQC_IBDQC_${Sample}_NoNA

#Remove Variant in R2 single
plink --bfile ${ana_dir}/temp/5.raw_snpQC_IBDQC_${Sample}_NoNA \
  --show-tags ${ana_dir}/temp/3.raw_snpQC_indepSNP.prune.in \
  --list-all \
 --tag-kb 250 \
 --tag-r2 0.8 \
 --threads ${Thread} \
 --allow-no-sex \
 --out ${ana_dir}/temp/6.raw_tagQC

awk '{print $1, $8}' ${ana_dir}/temp/6.raw_tagQC.tags.list | grep 'NONE' | awk '{print $1}' > ${ana_dir}/temp/7.raw_exclude.list

plink --bfile ${ana_dir}/temp/5.raw_snpQC_IBDQC_${Sample}_NoNA \
 --exclude ${ana_dir}/temp/7.raw_exclude.list \
 --make-bed \
 --threads ${Thread} \
 --allow-no-sex \
 --out ${ana_dir}/temp/8.raw_${Sample}_NoNA

#prune for MDS and PCA
plink --bfile ${ana_dir}/temp/8.raw_${Sample}_NoNA \
 --threads ${Thread} \
 --exclude ${inversion} \
 --range \
 --indep-pairwise 50 5 0.2 \
 --allow-no-sex \
 --out ${ana_dir}/temp/9.raw_indepSNP

plink --bfile ${ana_dir}/temp/8.raw_${Sample}_NoNA \
 --extract ${ana_dir}/temp/9.raw_indepSNP.prune.in \
 --threads ${Thread} \
 --genome \
 --allow-no-sex \
 --out ${ana_dir}/temp/10.raw_${Sample}_NoNA_genome

#MDS
mkdir ${ana_dir}/MDS
plink --bfile ${ana_dir}/temp/8.raw_${Sample}_NoNA \
 --extract ${ana_dir}/temp/9.raw_indepSNP.prune.in \
 --threads ${Thread} \
 --make-bed \
 --allow-no-sex \
 --out ${ana_dir}/temp/11.raw_indep_${Sample}_NoNA


