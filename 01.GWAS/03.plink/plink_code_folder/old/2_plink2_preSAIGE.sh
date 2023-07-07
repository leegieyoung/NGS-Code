##!/bin/sh
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
Case_pheno="${ana_dir}/temp/11.Case_pheno.txt"
Control_pheno="${ana_dir}/temp/11.Control_pheno.txt"
inversion="${Dir}/REFERENCE/inversion.txt"
Rcode="/mnt/nas/gylee/0.GWAS/Code/R"
p2_Dir="/mnt/nas/gylee/Singurality/plink2"
##=================================================================

mkdir -p ${Dir}/2.plink_result/NoImputed_${Sample}/temp
mkdir -p ${Dir}/2.plink_result/NoImputed_${Sample}/result
ana_dir="${Dir}/2.plink_result/NoImputed_${Sample}"

##=========What is Code ? =======================================
echo "${code_path}/2_Single_version" > ${ana_dir}/2_single_version

${p2_Dir}/plink2 --bfile ${MaMi_dir}/QC_Imputed_${Sample} \
 --glm \
 --threads ${Thread} \
 --memory 240000 \
 --out ${ana_dir}/temp/1.raw_${Sample}_assoc

#grep -v "NA" ${ana_dir}/raw_${Sample}_assoc.assoc | grep -v "CHR" | awk '{print $2}' > ${ana_dir}/extract_SNP_list.txt
awk '$9 != "NA" {print $0}' ${ana_dir}/temp/1.raw_${Sample}_assoc.assoc | grep -v "CHR" | awk '{print $2}' > ${ana_dir}/extract_SNP_list.txt


#==========================================================

${p2_Dir}/plink2 --bfile ${MaMi_dir}/QC_Imputed_${Sample} \
 --extract ${ana_dir}/extract_SNP_list.txt \
 --make-bed \
 --memory 240000 \
 --threads ${Thread} \
 --out ${ana_dir}/temp/2.raw_snpQC_${Sample}_NoNA
#==========================================================

#prune
${p2_Dir}/plink2 --bfile ${ana_dir}/temp/2.raw_snpQC_${Sample}_NoNA \
 --exclude ${inversion} \
 --memory 240000 \
 --indep-pairwise 50 5 0.02 \
 --threads ${Thread} \
 --out ${ana_dir}/temp/3.raw_snpQC_indepSNP

${p2_Dir}/plink2 --bfile ${ana_dir}/temp/2.raw_snpQC_${Sample}_NoNA \
 --extract ${ana_dir}/temp/3.raw_snpQC_indepSNP.prune.in \
 --make-bed \
 --memory 240000 \
 --threads ${Thread} \
 --out ${ana_dir}/temp/8.raw_${Sample}_NoNA


