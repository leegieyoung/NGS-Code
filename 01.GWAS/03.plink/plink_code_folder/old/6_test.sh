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
Case_pheno="${ana_dir}/temp/11.Case_pheno.txt"
Control_pheno="${ana_dir}/temp/11.Control_pheno.txt"
inversion="${Dir}/REFERENCE/inversion.txt"
Rcode="/mnt/nas/gylee/0.GWAS/Code/R"
p2_Dir="/mnt/nas/gylee/Singurality/plink2"
#=================================================================

mkdir -p ${Dir}/2.plink_result/NoImputed_${Sample}/temp
mkdir -p ${Dir}/2.plink_result/NoImputed_${Sample}/result
ana_dir="${Dir}/2.plink_result/NoImputed_${Sample}"

#=========What is Code ? =======================================
echo "${code_path}/2_Single_version" > ${ana_dir}/2_single_version


#1. calculate allele counts for each marker in the large plink file with hard called genotypes
#${p2_Dir}/plink2 --bfile ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC_${Sample} \
# --freq counts \
# --threads ${Thread} \
# --out ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC_${Sample}_freq

#2 randomly extract IDs for markers falling in the two MAC categories
#cat <(tail -n +2 ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC_${Sample}_freq.acount | awk '(($6-$5) < 20 && ($6-$5) >= 10) || ($5 < 20 && $5 >= 10) {print $2}' | shuf -n 1000) \
# <(tail -n +2 ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC_${Sample}_freq.acount | awk ' $5 >= 20 && ($6-$5)>= 20 {print $2}' | shuf -n 1000) > ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC_${Sample}.markerid.list

#None
cat <(tail -n +2 ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC_${Sample}_freq.acount | awk '(($7-$6) < 20 && ($7-$6) >= 10) || ($6 < 20 && $6 >= 10) {print $2}' | shuf -n 1000) > ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC_${Sample}.markerid.list1
#1000
cat <(tail -n +2 ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC_${Sample}_freq.acount | awk ' $6 >= 20 && ($7-$6)>= 20 {print $2}' | shuf -n 1000) > ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC_${Sample}.markerid.list2

#3. extract markers from the large plink file
#${p2_Dir}/plink2 --bfile ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC_${Sample} \
# --extract ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC_${Sample}.markerid.list \
# --make-bed \
# --threads ${Thread} \
# --out ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/1000_for_vr_${Sample}
