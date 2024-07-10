#!/bin/sh
if [ $# -ne 1 ];then
        echo "Please enter Sample_Name"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
Sample=$1
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

cp -v ${beagle_folder}/imputed_QC_${Sample}/* ${Sample_folder}/

for Chrom in {1..22}
do
plink --vcf ${Sample_folder}/CHR${Chrom}_QC_${Sample}_noIndel_rmDupID_beagle.vcf.gz --make-bed --out ${result_folder}/raw_merge/NoQC_CHR${Chrom}_${Sample}
done

for Chrom in {2..22}
do
echo ${result_folder}/raw_merge/NoQC_CHR${Chrom}_${Sample}.bed ${result_folder}/raw_merge/NoQC_CHR${Chrom}_${Sample}.bim ${result_folder}/raw_merge/NoQC_CHR${Chrom}_${Sample}.fam > ${result_folder}/raw_merge/mergelist_${Chrom}.txt
done
cat ${result_folder}/raw_merge/mergelist_*.txt > ${result_folder}/raw_merge/mergelist.txt

plink --bfile ${result_folder}/raw_merge/NoQC_CHR1_${Sample} --merge-list ${result_folder}/raw_merge/mergelist.txt --make-bed --out ${result_folder}/QC/raw_NoQC_Imputed_${Sample} 

# QC
plink --bfile ${result_folder}/QC/raw_NoQC_Imputed_${Sample} \
 --geno 0.05 \
 --make-bed \
 --out ${result_folder}/QC/raw_QC_Imputed_${Sample}_g

plink --bfile ${result_folder}/QC/raw_QC_Imputed_${Sample}_g \
 --mind 0.05 \
 --make-bed \
 --out ${result_folder}/QC/raw_QC_Imputed_${Sample}_g_m

plink --bfile ${result_folder}/QC/raw_QC_Imputed_${Sample}_g_m \
 --maf 0.001 \
 --make-bed \
 --out ${result_folder}/QC/raw_QC_Imputed_${Sample}_g_m_maf

plink --bfile ${result_folder}/QC/raw_QC_Imputed_${Sample}_g_m_maf \
 --hwe 1e-06 \
 --make-bed \
 --out ${result_folder}/QC/raw_QC_Imputed_${Sample}_g_m_maf_hwe

cp ${result_folder}/QC/raw_QC_Imputed_${Sample}_g_m_maf_hwe.bed ${result_folder}/MaMi/QC_Imputed_${Sample}.bed
cp ${result_folder}/QC/raw_QC_Imputed_${Sample}_g_m_maf_hwe.bim ${result_folder}/MaMi/QC_Imputed_${Sample}.bim
cp /scratch/x1997a11/GWAS/pdxen_AD/Sample_folder/${Sample}/${Sample}.fam  ${result_folder}/MaMi/QC_Imputed_${Sample}.fam

#rsID change
echo "=========================================="
echo "
          Change rsID that CHR:POSI"
echo "=============================+============"
sh ${code_path}/3_result_MaMi_rsIDchange_chrposi.sh Imputed_${Sample}
