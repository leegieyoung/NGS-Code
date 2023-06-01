#!/bin/sh
if [ $# -ne 2 ];then
        echo "Please enter Sample_Name"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
Sample=$1
Chrom=$2

GWAS_path="/scratch/x1997a11/GWAS/pdxen_AD/result_folder"
code_path="/scratch/x1997a11/GWAS/pdxen_AD/Code_folder"
mkdir ${GWAS_path}/QC_Imputed_${Sample}
result_folder="${GWAS_path}/QC_Imputed_${Sample}"
#Case 원본파일이 담긴 곳
mkdir /scratch/x1997a11/GWAS/pdxen_AD/Sample_folder/Imputed_${Sample}
mkdir /scratch/x1997a11/GWAS/pdxen_AD/Sample_folder/Imputed_${Sample}/QC_Imputed_${Sample}
Sample_folder="/scratch/x1997a11/GWAS/pdxen_AD/Sample_folder/Imputed_${Sample}/QC_Imputed_${Sample}"
beagle_folder="/scratch/x1997a11/GWAS/pdxen_AD/beagle/Result_folder"

#1.SampleQC(geno mind impute-sex hwe)

mkdir ${result_folder}/raw_merge
mkdir ${result_folder}/QC
mkdir ${result_folder}/MaMi

#
echo ''${Chrom}_'vcf to '${Chrom}_'bfile'
cp -v ${beagle_folder}/imputed_${Sample}/CHR${Chrom}_${Sample}_noIndel_rmDupID_beagle.vcf.gz ${Sample_folder}/CHR${Chrom}_${Sample}_noIndel_rmDupID_beagle.vcf.gz
plink --vcf ${Sample_folder}/CHR${Chrom}_${Sample}_noIndel_rmDupID_beagle.vcf.gz  --make-bed --keep-allele-order --double-id --out ${result_folder}/raw_merge/NoQC_CHR${Chrom}_${Sample}

#rmIDdot
echo ''${Chrom}_'rmIDdot'
awk '$2 != "." && $5 != "." && $6 != "." {print $2}' ${result_folder}/raw_merge/NoQC_CHR${Chrom}_${Sample}.bim > ${result_folder}/raw_merge/NoQC_rmIDdot_CHR${Chrom}_${Sample}.txt
plink --bfile ${result_folder}/raw_merge/NoQC_CHR${Chrom}_${Sample} --extract ${result_folder}/raw_merge/NoQC_rmIDdot_CHR${Chrom}_${Sample}.txt --make-bed --keep-allele-order --out ${result_folder}/raw_merge/NoQC_rmIDdot_CHR${Chrom}_${Sample}
