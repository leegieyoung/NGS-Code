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
plink --bfile ${Sample_folder}/${Sample} --recode vcf-iid --out ${result_folder}/QC/raw_Case
cat ${result_folder}/QC/raw_Case.vcf | awk '$4 != "-" && $5 != "-" && $4 != "I" && $5 != "I" && $4 != "D" && $5 != "D" && $4 != "0" && $5 != "0" && $4 != "." && $5 != "."' > ${result_folder}/QC/raw_Case_noIndel.vcf
bcftools norm -d all -O v -o ${result_folder}/QC/XYMT_raw_Case_noIndel_rmDupID.vcf ${result_folder}/QC/raw_Case_noIndel.vcf
grep -v -w '^23' ${result_folder}/QC/XYMT_raw_Case_noIndel_rmDupID.vcf | grep -v -w '^24' | grep -v -w '^25' | grep -v -w '^26' > ${result_folder}/QC/raw_Case_noIndel_rmDupID.vcf

plink --vcf ${result_folder}/QC/raw_Case_noIndel_rmDupID.vcf --make-bed --out ${result_folder}/QC/raw_Case

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


echo "=========================================="
echo "             

                    Sample QC_ing... 

"
echo "=========================================+"
#impute-sex 를 사용하지 않는 경우
plink --bfile ${result_folder}/QC/raw_Case --geno 0.2 --make-bed --out ${result_folder}/QC/raw_Case_g
plink --bfile ${result_folder}/QC/raw_Case_g --mind 0.2 --make-bed --out ${result_folder}/QC/raw_Case_g_m
#plink --bfile ${result_folder}/QC/raw_Case_g_m --impute-sex --make-bed --out ${result_folder}/QC/raw_Case_g_m_sex

#cp ${result_folder}/QC/raw_Case_g_m_sex.fam ${result_folder}/QC/SEXerror_raw_Case_g_m_sex.fam
#awk '{sub(/0/,"2",$5);print > "'${result_folder}'/QC/raw_Case_g_m_sex.fam"}' ${result_folder}/QC/SEXerror_raw_Case_g_m_sex.fam
plink --bfile ${result_folder}/QC/raw_Case_g_m \
 --maf 0.01 \
 --make-bed \
 --out ${result_folder}/QC/raw_Case_g_m_maf

plink --bfile ${result_folder}/QC/raw_Case_g_m_maf \
 --hwe 1e-6 \
 --make-bed \
 --out ${result_folder}/QC/raw_Case_g_m_maf_hwe

echo "=========================================="
echo "             

                    Sample QC_end 

"
echo "=========================================+"


#===============================================================================================================================

#Create QC & MaMi/MiMa.binary file 

cp -v ${result_folder}/QC/raw_Case_g_m_maf_hwe.bim ${result_folder}/MaMi/QC_${Sample}.bim
cp -v ${result_folder}/QC/raw_Case_g_m_maf_hwe.bed ${result_folder}/MaMi/QC_${Sample}.bed
awk '{sub(/-9/,"2",$6);print > "'${result_folder}'/MaMi/pheCase_QC_'${Sample}.fam'"}' ${result_folder}/QC/raw_Case_g_m_maf_hwe.fam
awk '{sub(/0/,"1",$5);print > "'${result_folder}'/MaMi/QC_'${Sample}.fam'"}' ${result_folder}/MaMi/pheCase_QC_${Sample}.fam
#paste -d '\t' ${result_folder}/MiMa/raw_12345bim.bim  ${result_folder}/QC/raw_Minor_merge_case.bim > ${result_folder}/MiMa/QC_${Sample}.bim
#cp -v ${result_folder}/QC/raw_Case_g_m_maf_hwe.bed ${result_folder}/MiMa/QC_${Sample}.bed
#cp -v ${result_folder}/QC/raw_Case_g_m_maf_hwe.fam ${result_folder}/MiMa/QC_${Sample}.fam
#=============================================================================
#분석에 사용한 데이터 제거
#rm ${result_folder}/QC/raw_*
#rm ${result_folder}/MaMi/raw_*
#rm ${result_folder}/MiMa/raw_*

