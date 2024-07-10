#!/bin/sh

sample_path="/scratch/x1997a11/GWAS/pdxen_AD/Sample_folder/Imputed_AD_199/QC_Imputed_AD_199/AD_14_AD_185bfile"
keep_path="/scratch/x1997a11/GWAS/pdxen_AD/Sample_folder/repeat_check"
result_path="/scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check"
mkdir ${result_path}/Case
mkdir ${result_path}/Control
mkdir ${result_path}/Case_Control
#Case
for A in $(seq 1 2)
do
mkdir ${result_path}/Case/Case$A
mkdir ${result_path}/Control/Control$A

paste ${keep_path}/Case$A.txt -d '\t' ${keep_path}/Case$A.txt > ${result_path}/Case/Case$A/Case$A.txt
sed -i 's///g' ${result_path}/Case/Case$A/Case$A.txt
paste ${keep_path}/Control$A.txt -d '\t' ${keep_path}/Control$A.txt > ${result_path}/Control/Control$A/Control$A.txt
sed -i 's///g' ${result_path}/Control/Control$A/Control$A.txt

sh logi_repeat_check.sh $A 
done




