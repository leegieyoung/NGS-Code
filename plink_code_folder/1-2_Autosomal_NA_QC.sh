#!/bin/sh
if [ $# -ne 1 ];then
        echo "Please enter Sample_Name"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
#beagle 진행시Same position 제거를 안했을 경우 사용하는 코드
Sample=$1
Dir="/mnt/nas/gylee/0.GWAS"
impute_dir="${Dir}/1.Input/${Sample}/0.Impute"
QC_dir="${Dir}/1.Input/${Sample}/1.QC"
MaMi_dir="${Dir}/1.Input/${Sample}/2.MaMi"

#Remove Unknwon variant 
awk '$1 != "26" && $1 != "25" && $1 != "24" && $1 != "23" {print $0}' ${QC_dir}/${Sample}_g_m_maf_hwe.bim | awk '$5 != "D" && $6 != "D" && $5 != "I" && $6 != "I" && !index($5, "<") && !index($6, "<") {print $2}' > ${QC_dir}/${Sample}_QC_auto_misVari.txt

plink --bfile ${QC_dir}/${Sample}_g_m_maf_hwe \
	--extract ${QC_dir}/${Sample}_QC_auto_misVari.txt \
	--keep-allele-order \
	--allow-no-sex \
	--make-bed \
	--out ${MaMi_dir}/QC_Imputed_${Sample}


echo ""
echo "===================Delete Used Data========================"
echo ""
