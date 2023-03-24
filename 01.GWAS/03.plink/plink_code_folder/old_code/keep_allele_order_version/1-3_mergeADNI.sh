#!/bin/sh

GWAS_path="/scratch/x1997a11/GWAS/pdxen_AD/result_folder"
Code_path="/scratch/x1997a11/GWAS/pdxen_AD/Code_folder/keep_allele_order_version"
mkdir ${GWAS_path}/QC_Imputed_mergeADNI
mkdir ${GWAS_path}/QC_Imputed_mergeADNI/MaMi
mkdir ${GWAS_path}/QC_Imputed_mergeADNI/QC

echo ${GWAS_path}/QC_Imputed_ADNI_GO2/MaMi/QC_Imputed_ADNI_GO2.bed ${GWAS_path}/QC_Imputed_ADNI_GO2/MaMi/QC_Imputed_ADNI_GO2.bim ${GWAS_path}/QC_Imputed_ADNI_GO2/MaMi/QC_Imputed_ADNI_GO2.fam > ${GWAS_path}/QC_Imputed_ADNI_GO2/MaMi/QC_Imputed_ADNI_GO2.mergelist

echo ${GWAS_path}/QC_Imputed_ADNI3/MaMi/QC_Imputed_ADNI3.bed ${GWAS_path}/QC_Imputed_ADNI3/MaMi/QC_Imputed_ADNI3.bim ${GWAS_path}/QC_Imputed_ADNI3/MaMi/QC_Imputed_ADNI3.fam > ${GWAS_path}/QC_Imputed_ADNI3/MaMi/QC_Imputed_ADNI3.mergelist

cat ${GWAS_path}/QC_Imputed_ADNI_GO2/MaMi/QC_Imputed_ADNI_GO2.mergelist ${GWAS_path}/QC_Imputed_ADNI3/MaMi/QC_Imputed_ADNI3.mergelist > ${GWAS_path}/QC_Imputed_mergeADNI/MaMi/mergelist.txt

plink --bfile /scratch/x1997a11/GWAS/pdxen_AD/result_folder/QC_Imputed_ADNI1/MaMi/QC_Imputed_ADNI1 \
 --merge-list ${GWAS_path}/QC_Imputed_mergeADNI/MaMi/mergelist.txt \
 --make-bed \
 --keep-allele-order \
 --allow-no-sex \
 --out ${GWAS_path}/QC_Imputed_mergeADNI/MaMi/QC_Imputed_mergeADNI

#missnp 제거
plink --bfile ${GWAS_path}/QC_Imputed_ADNI1/MaMi/QC_Imputed_ADNI1 \
 --exclude ${GWAS_path}/QC_Imputed_mergeADNI/MaMi/QC_Imputed_mergeADNI-merge.missnp \
 --keep-allele-order \
 --make-bed \
 --allow-no-sex \
 --out ${GWAS_path}/QC_Imputed_mergeADNI/MaMi/QC_Imputed_rmmissnp_ADNI1

plink --bfile ${GWAS_path}/QC_Imputed_ADNI_GO2/MaMi/QC_Imputed_ADNI_GO2 \
 --exclude ${GWAS_path}/QC_Imputed_mergeADNI/MaMi/QC_Imputed_mergeADNI-merge.missnp \
 --keep-allele-order \
 --make-bed \
 --allow-no-sex \
 --out ${GWAS_path}/QC_Imputed_mergeADNI/MaMi/QC_Imputed_rmmissnp_ADNI_GO2

plink --bfile ${GWAS_path}/QC_Imputed_ADNI3/MaMi/QC_Imputed_ADNI3 \
 --exclude ${GWAS_path}/QC_Imputed_mergeADNI/MaMi/QC_Imputed_mergeADNI-merge.missnp \
 --keep-allele-order \
 --make-bed \
 --allow-no-sex \
 --out ${GWAS_path}/QC_Imputed_mergeADNI/MaMi/QC_Imputed_rmmissnp_ADNI3

echo ${GWAS_path}/QC_Imputed_mergeADNI/MaMi/QC_Imputed_rmmissnp_ADNI_GO2.bed ${GWAS_path}/QC_Imputed_mergeADNI/MaMi/QC_Imputed_rmmissnp_ADNI_GO2.bim ${GWAS_path}/QC_Imputed_mergeADNI/MaMi/QC_Imputed_rmmissnp_ADNI_GO2.fam > ${GWAS_path}/QC_Imputed_mergeADNI/MaMi/QC_Imputed_rmmissnp_ADNI_GO2.mergelist

echo ${GWAS_path}/QC_Imputed_mergeADNI/MaMi/QC_Imputed_rmmissnp_ADNI3.bed ${GWAS_path}/QC_Imputed_mergeADNI/MaMi/QC_Imputed_rmmissnp_ADNI3.bim ${GWAS_path}/QC_Imputed_mergeADNI/MaMi/QC_Imputed_rmmissnp_ADNI3.fam > ${GWAS_path}/QC_Imputed_mergeADNI/MaMi/QC_Imputed_rmmissnp_ADNI3.mergelist

cat ${GWAS_path}/QC_Imputed_mergeADNI/MaMi/QC_Imputed_rmmissnp_ADNI_GO2.mergelist ${GWAS_path}/QC_Imputed_mergeADNI/MaMi/QC_Imputed_rmmissnp_ADNI3.mergelist > ${GWAS_path}/QC_Imputed_mergeADNI/MaMi/rmmissnp.mergelist

#merge
plink --bfile ${GWAS_path}/QC_Imputed_mergeADNI/MaMi/QC_Imputed_rmmissnp_ADNI1 \
 --merge-list ${GWAS_path}/QC_Imputed_mergeADNI/MaMi/rmmissnp.mergelist \
 --keep-allele-order \
 --make-bed \
 --allow-no-sex \
 --out ${GWAS_path}/QC_Imputed_mergeADNI/QC/NoQC_Imputed_rmmissnp_mergeADNI

plink --bfile ${GWAS_path}/QC_Imputed_mergeADNI/QC/NoQC_Imputed_rmmissnp_mergeADNI \
 --keep-allele-order \
 --make-bed \
 --allow-no-sex \
 --geno 0.2 \
 --out ${GWAS_path}/QC_Imputed_mergeADNI/QC/NoQC_Imputed_rmmissnp_mergeADNI_g 

plink --bfile ${GWAS_path}/QC_Imputed_mergeADNI/QC/NoQC_Imputed_rmmissnp_mergeADNI_g \
 --keep-allele-order \
 --make-bed \
 --allow-no-sex \
 --mind 0.2 \
 --out ${GWAS_path}/QC_Imputed_mergeADNI/QC/NoQC_Imputed_rmmissnp_mergeADNI_g_m

plink --bfile ${GWAS_path}/QC_Imputed_mergeADNI/QC/NoQC_Imputed_rmmissnp_mergeADNI_g_m \
 --keep-allele-order \
 --make-bed \
 --allow-no-sex \
 --maf 0.01 \
 --out ${GWAS_path}/QC_Imputed_mergeADNI/QC/NoQC_Imputed_rmmissnp_mergeADNI_g_m_maf 

plink --bfile ${GWAS_path}/QC_Imputed_mergeADNI/QC/NoQC_Imputed_rmmissnp_mergeADNI_g_m_maf \
 --keep-allele-order \
 --make-bed \
 --allow-no-sex \
 --hwe 1e-6 \
 --out ${GWAS_path}/QC_Imputed_mergeADNI/QC/NoQC_Imputed_rmmissnp_mergeADNI_g_m_maf_hwe

cp ${GWAS_path}/QC_Imputed_mergeADNI/QC/NoQC_Imputed_rmmissnp_mergeADNI_g_m_maf_hwe.bed ${GWAS_path}/QC_Imputed_mergeADNI/MaMi/QC_Imputed_rmmissnp_mergeADNI.bed
cp ${GWAS_path}/QC_Imputed_mergeADNI/QC/NoQC_Imputed_rmmissnp_mergeADNI_g_m_maf_hwe.bim ${GWAS_path}/QC_Imputed_mergeADNI/MaMi/QC_Imputed_rmmissnp_mergeADNI.bim
cp ${GWAS_path}/QC_Imputed_mergeADNI/QC/NoQC_Imputed_rmmissnp_mergeADNI_g_m_maf_hwe.fam ${GWAS_path}/QC_Imputed_mergeADNI/MaMi/QC_Imputed_rmmissnp_mergeADNI.fam
