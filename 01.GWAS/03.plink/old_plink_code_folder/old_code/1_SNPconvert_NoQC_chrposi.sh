#!/bin/sh
if [ $# -ne 1 ];then
        echo "Please enter Sample_Name"
               exit
fi
Sample=$1

GWAS_path="/scratch/x1997a11/GWAS/pdxen_AD/result_folder"
code_path="/scratch/x1997a11/GWAS/pdxen_AD/Code_folder"
mkdir ${GWAS_path}/QC_${Sample}
result_folder="${GWAS_path}/QC_${Sample}"
#Case 원본파일이 담긴 곳
Sample_folder="/scratch/x1997a11/GWAS/pdxen_AD/Sample_folder/${Sample}"
pheno="/scratch/x1997a11/GWAS/pdxen_AD/reference_folder/${Sample}.txt"


#1.SampleQC(geno mind impute-sex hwe)

#Case data 작업폴더
mkdir ${result_folder}/QC
#=======================================================================

#Indel 제거
#plink --bfile ${Sample_folder}/${Sample} --recode vcf-iid --out ${result_folder}/QC/raw_Case
#cat ${result_folder}/QC/raw_Case.vcf | awk '$4 != "-" && $5 != "-" && $4 != "I" && $5 != "I" && $4 != "D" && $5 != "D" && $4 != "0" && $5 != "0" && $4 != "." && $5 != "."' > ${result_folder}/QC/raw_Case_noIndel.vcf
#bcftools norm -d all -O v -o ${result_folder}/QC/raw_Case_noIndel_rmDupID.vcf ${result_folder}/QC/raw_Case_noIndel.vcf

#plink --vcf ${result_folder}/QC/raw_Case_noIndel_rmDupID.vcf --make-bed --out ${result_folder}/QC/raw_Case

#X, Y, MT, XM 제거
cp ${Sample_folder}/${Sample}.bed ${result_folder}/QC/XYMT_raw_Case.bed
cp ${Sample_folder}/${Sample}.bim ${result_folder}/QC/XYMT_raw_Case.bim
cp ${Sample_folder}/${Sample}.fam ${result_folder}/QC/XYMT_raw_Case.fam

plink -bfile ${result_folder}/QC/XYMT_raw_Case --recode vcf-iid --out ${result_folder}/QC/XYMT_raw_Case
grep -v -w '^23' ${result_folder}/QC/XYMT_raw_Case.vcf | grep -v -w '^24' | grep -v -w '^25' | grep -v -w '^26' > ${result_folder}/QC/raw_Case.vcf
plink --vcf ${result_folder}/QC/raw_Case.vcf --make-bed --out ${result_folder}/QC/raw_Case

#chr:posi_Major/Minor SNP
awk '{$1=$2="";print $0}' ${result_folder}/QC/raw_Case.bim > ${result_folder}/QC/raw_other_col_merge_case.bim
awk '{print $1}' ${result_folder}/QC/raw_Case.bim > ${result_folder}/QC/raw_col_chr_merge_case.bim
awk '{print $4}' ${result_folder}/QC/raw_Case.bim > ${result_folder}/QC/raw_col_posi_merge_case.bim
awk '{print $6}' ${result_folder}/QC/raw_Case.bim > ${result_folder}/QC/raw_Major_merge_case.bim
awk '{print $5}' ${result_folder}/QC/raw_Case.bim > ${result_folder}/QC/raw_Minor_merge_case.bim

#chr:posi
paste -d : ${result_folder}/QC/raw_col_chr_merge_case.bim ${result_folder}/QC/raw_col_posi_merge_case.bim > ${result_folder}/QC/raw_12_case.bim

#make bim
mkdir ${result_folder}/MaMi
paste -d '\t' ${result_folder}/QC/raw_col_chr_merge_case.bim ${result_folder}/QC/raw_12_case.bim  > ${result_folder}/QC/raw_12bim.bim
awk '{print($0"\t'0'")}' ${result_folder}/QC/raw_12bim.bim > ${result_folder}/QC/raw_123bim.bim
paste -d '\t' ${result_folder}/QC/raw_123bim.bim ${result_folder}/QC/raw_col_posi_merge_case.bim > ${result_folder}/QC/raw_1234bim.bim
paste -d '\t' ${result_folder}/QC/raw_1234bim.bim ${result_folder}/QC/raw_Minor_merge_case.bim > ${result_folder}/QC/raw_12345bim.bim
paste -d '\t' ${result_folder}/QC/raw_12345bim.bim  ${result_folder}/QC/raw_Major_merge_case.bim > ${result_folder}/QC/raw_Case.bim

#===============================================================================================================================

#Create QC & MaMi/MiMa.binary file 

cp -v ${result_folder}/QC/raw_Case.bim ${result_folder}/MaMi/QC_${Sample}.bim
cp -v ${result_folder}/QC/raw_Case.bed ${result_folder}/MaMi/QC_${Sample}.bed
cp -v ${Sample_folder}/${Sample}.fam ${result_folder}/MaMi/QC_${Sample}.fam

#paste -d '\t' ${result_folder}/MiMa/raw_12345bim.bim  ${result_folder}/QC/raw_Minor_merge_case.bim > ${result_folder}/MiMa/QC_${Sample}.bim
#cp -v ${result_folder}/QC/raw_Case_g_m_hwe.bed ${result_folder}/MiMa/QC_${Sample}.bed
#cp -v ${result_folder}/QC/raw_Case_g_m_hwe.fam ${result_folder}/MiMa/QC_${Sample}.fam
#=============================================================================
#분석에 사용한 데이터 제거
#rm ${result_folder}/QC/raw_*
#rm ${result_folder}/MaMi/raw_*
#rm ${result_folder}/MiMa/raw_*

