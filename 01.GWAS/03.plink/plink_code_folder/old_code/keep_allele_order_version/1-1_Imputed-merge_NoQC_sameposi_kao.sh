#!/bin/sh
if [ $# -ne 1 ];then
        echo "Please enter Sample_Name"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
#beagle 진행시Same position 제거를 안했을 경우 사용하는 코드
Sample=$1
GWAS_path="/data/keeyoung/GWAS/result"
code_path="/scratch/x1997a11/GWAS/pdxen_AD/Code_folder"
mkdir ${GWAS_path}/QC_Imputed_${Sample}
result_folder="${GWAS_path}/QC_Imputed_${Sample}"
#imputed_sample folder
beagle_folder="/data/keeyoung/GWAS/beagle_result"
Sample_folder="${beagle_folder}/3.${Sample}_imputed_${Sample}"
raw_folder="/data/keeyoung/GWAS/raw_data/"
#beagle_folder="/scratch/x1997a11/GWAS/pdxen_AD/beagle/Result_folder"

#1.SampleQC(geno mind impute-sex hwe)
mkdir ${result_folder}/raw_merge
mkdir ${result_folder}/QC
mkdir ${result_folder}/MaMi

cp -v ${beagle_folder}/imputed_QC_${Sample}/* ${Sample_folder}/
for Chrom in {1..22}
do
plink --vcf ${Sample_folder}/CHR${Chrom}_${Sample}_noIndel_rmDupID_beagle.vcf.gz --double-id --make-bed --keep-allele-order --out ${result_folder}/raw_merge/NoQC_CHR${Chrom}_${Sample}
done

rmIDdot
for Chrom in {1..22}
do
echo ''${Chrom}_'rmIDdot'
awk '$2 != "." && $5 != "." && $6 != "." {print $2}' ${result_folder}/raw_merge/NoQC_CHR${Chrom}_${Sample}.bim > ${result_folder}/raw_merge/NoQC_rmIDdot_CHR${Chrom}_${Sample}.txt
plink --bfile ${result_folder}/raw_merge/NoQC_CHR${Chrom}_${Sample} --extract ${result_folder}/raw_merge/NoQC_rmIDdot_CHR${Chrom}_${Sample}.txt --make-bed --keep-allele-order --out ${result_folder}/raw_merge/NoQC_rmIDdot_CHR${Chrom}_${Sample}
done

for Chrom in {2..22}
do
echo ${result_folder}/raw_merge/NoQC_rmIDdot_CHR${Chrom}_${Sample}.bed ${result_folder}/raw_merge/NoQC_rmIDdot_CHR${Chrom}_${Sample}.bim ${result_folder}/raw_merge/NoQC_rmIDdot_CHR${Chrom}_${Sample}.fam > ${result_folder}/raw_merge/mergelist_${Chrom}.txt
done
cat ${result_folder}/raw_merge/mergelist_*.txt > ${result_folder}/raw_merge/mergelist.txt

plink --bfile ${result_folder}/raw_merge/NoQC_rmIDdot_CHR1_${Sample} \
 --merge-list ${result_folder}/raw_merge/mergelist.txt \
 --make-bed \
 --keep-allele-order \
 --out ${result_folder}/QC/raw_NoQC_same_posi_Imputed_${Sample} 


#rm Same position
grep '^Warning' ${result_folder}/QC/raw_NoQC_same_posi_Imputed_${Sample}.log > ${result_folder}/QC/same_position.log
sed -i "s/'/\t/g" ${result_folder}/QC/same_position.log
awk '{print $3}' ${result_folder}/QC/same_position.log > ${result_folder}/QC/same_position.rsid
plink --bfile ${result_folder}/QC/raw_NoQC_same_posi_Imputed_${Sample} \
 --exclude ${result_folder}/QC/same_position.rsid \
 --keep-allele-order \
 --make-bed \
 --out ${result_folder}/QC/raw_NoQC_Imputed_${Sample}
