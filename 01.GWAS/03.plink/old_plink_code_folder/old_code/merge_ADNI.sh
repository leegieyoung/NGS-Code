#!/bin/sh
path="/scratch/x1997a11/GWAS/pdxen_AD/Code_folder"
GWAS_path="/scratch/x1997a11/GWAS/pdxen_AD/result_folder"

echo ${GWAS_path}/QC_Imputed_ADNI_GO2/MaMi/QC_Imputed_rmmissnp_ADNI_GO2.bed ${GWAS_path}/QC_Imputed_ADNI_GO2/MaMi/QC_Imputed_rmmissnp_ADNI_GO2.bim ${GWAS_path}/QC_Imputed_ADNI_GO2/MaMi/QC_Imputed_rmmissnp_ADNI_GO2.fam > ${GWAS_path}/QC_Imputed_ADNI_GO2/MaMi/mergelist_ADNI_GO2.txt
echo ${GWAS_path}/QC_Imputed_ADNI3/MaMi/QC_Imputed_rmmissnp_ADNI3.bed ${GWAS_path}/QC_Imputed_ADNI3/MaMi/QC_Imputed_rmmissnp_ADNI3.bim ${GWAS_path}/QC_Imputed_ADNI3/MaMi/QC_Imputed_rmmissnp_ADNI3.fam > ${GWAS_path}/QC_Imputed_ADNI3/MaMi/mergelist_ADNI3.txt
cat ${GWAS_path}/QC_Imputed_ADNI_GO2/MaMi/mergelist_ADNI_GO2.txt ${GWAS_path}/QC_Imputed_ADNI3/MaMi/mergelist_ADNI3.txt > ${GWAS_path}/QC_Imputed_ADNI1/MaMi/mergelist.txt
plink --bfile ${GWAS_path}/QC_Imputed_ADNI1/MaMi/QC_Imputed_rmmissnp_ADNI1 \
 --merge-list ${GWAS_path}/QC_Imputed_ADNI1/MaMi/mergelist.txt \
 --make-bed \
 --keep-allele-order \
 --out ${GWAS_path}/QC_Imputed_ADNI1/MaMi/QC_Imputed_rmmissnp_mergeADNI
