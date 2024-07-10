#!/bin/sh
Sample_path="/scratch/x1997a11/GWAS/pdxen_AD/result_folder/Imputed_AD_14_Imputed_AD_185_analysis_folder/merge"
Result_folder="/scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check/plink_pphe"

Start=$1
num=99
End=$((Start+num))


for A in $(seq $Start $End)
do
B=$((A+2))
mkdir ${Result_folder}/check$A
awk '{print $1, $2, $3, $4, $5}' ${Sample_path}/Imputed_AD_14_Imputed_AD_185_NoNA.fam > ${Result_folder}/check$A/Nophe_check$A.fam
awk '{print $'${B}'}' ${Sample_path}/Imputed_AD_14_Imputed_AD_185_NoNA_pphe.pphe > ${Result_folder}/check$A/Phe_check$A.txt
paste ${Result_folder}/check$A/Nophe_check$A.fam ${Result_folder}/check$A/Phe_check$A.txt > ${Result_folder}/check$A/check$A.fam
plink \
 --bed ${Sample_path}/Imputed_AD_14_Imputed_AD_185_NoNA.bed \
 --bim ${Sample_path}/Imputed_AD_14_Imputed_AD_185_NoNA.bim \
 --fam ${Result_folder}/check$A/check$A.fam \
 --ci 0.95 \
 --covar ${Sample_path}/logistic.covar_mds.txt \
 --logistic \
 --hide-covar \
 --out ${Result_folder}/check$A/check${A}
done


#for A in $(seq 1 10)
#do
#B=$((A+2))
#echo $B
#done
