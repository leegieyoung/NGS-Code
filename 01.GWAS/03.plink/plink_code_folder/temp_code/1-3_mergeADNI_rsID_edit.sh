#!/bin/sh

output=$1
Sample1=$2
Sample2=$3

GWAS_path="/scratch/x1997a11/GWAS/pdxen_AD/result_folder"
Code_path="/scratch/x1997a11/GWAS/pdxen_AD/Code_folder/keep_allele_order_version"
mkdir ${GWAS_path}/QC_${output}
mkdir ${GWAS_path}/QC_${output}/MaMi
mkdir ${GWAS_path}/QC_${output}/QC

echo ${GWAS_path}/QC_${Sample2}/MaMi/QC_${Sample2}.bed ${GWAS_path}/QC_${Sample2}/MaMi/QC_${Sample2}.bim ${GWAS_path}/QC_${Sample2}/MaMi/QC_${Sample2}.fam > ${GWAS_path}/QC_${Sample2}/MaMi/QC_${Sample2}.mergelist


#echo ${GWAS_path}/QC_ADNI3/MaMi/QC_ADNI3.bed ${GWAS_path}/QC_ADNI3/MaMi/QC_ADNI3.bim ${GWAS_path}/QC_ADNI3/MaMi/QC_ADNI3.fam > ${GWAS_path}/QC_ADNI3/MaMi/QC_ADNI3.mergelist

cat ${GWAS_path}/QC_${Sample2}/MaMi/QC_${Sample2}.mergelist > ${GWAS_path}/QC_${output}/MaMi/mergelist.txt

plink --bfile ${GWAS_path}/QC_${Sample1}/MaMi/QC_${Sample1} \
 --merge-list ${GWAS_path}/QC_${output}/MaMi/mergelist.txt \
 --make-bed \
 --keep-allele-order \
 --allow-no-sex \
 --out ${GWAS_path}/QC_${output}/MaMi/QC_${output}

awk '{print $2}' ${GWAS_path}/QC_${output}/MaMi/QC_${output}.fam > ${GWAS_path}/QC_${output}/MaMi/mv_QC_${output}-merge.missnp
mv ${GWAS_path}/QC_${output}/MaMi/QC_${output}-merge.missnp ${GWAS_path}/QC_${output}/MaMi/mv_QC_${output}-merge.missnp


#missnp 제거
plink --bfile ${GWAS_path}/QC_${Sample1}/MaMi/QC_${Sample1} \
 --exclude ${GWAS_path}/QC_${output}/MaMi/mv_QC_${output}-merge.missnp \
 --keep-allele-order \
 --make-bed \
 --allow-no-sex \
 --out ${GWAS_path}/QC_${output}/MaMi/QC_rmmissnp_${Sample1}

plink --bfile ${GWAS_path}/QC_${Sample2}/MaMi/QC_${Sample2} \
 --exclude ${GWAS_path}/QC_${output}/MaMi/mv_QC_${output}-merge.missnp \
 --keep-allele-order \
 --make-bed \
 --allow-no-sex \
 --out ${GWAS_path}/QC_${output}/MaMi/QC_rmmissnp_${Sample2}

echo ${GWAS_path}/QC_${output}/MaMi/QC_rmmissnp_${Sample2}.bed ${GWAS_path}/QC_${output}/MaMi/QC_rmmissnp_${Sample2}.bim ${GWAS_path}/QC_${output}/MaMi/QC_rmmissnp_${Sample2}.fam > ${GWAS_path}/QC_${output}/MaMi/QC_rmmissnp_${Sample2}.mergelist


cat ${GWAS_path}/QC_${output}/MaMi/QC_rmmissnp_${Sample2}.mergelist > ${GWAS_path}/QC_${output}/MaMi/rmmissnp.mergelist

#merge
plink --bfile ${GWAS_path}/QC_${output}/MaMi/QC_rmmissnp_${Sample1} \
 --merge-list ${GWAS_path}/QC_${output}/MaMi/rmmissnp.mergelist \
 --keep-allele-order \
 --make-bed \
 --allow-no-sex \
 --out ${GWAS_path}/QC_${output}/QC/NoQC_rmmissnp_${output}

plink --bfile ${GWAS_path}/QC_${output}/QC/NoQC_rmmissnp_${output} \
 --keep-allele-order \
 --make-bed \
 --allow-no-sex \
 --geno 0.2 \
 --out ${GWAS_path}/QC_${output}/QC/NoQC_rmmissnp_${output}_g 

plink --bfile ${GWAS_path}/QC_${output}/QC/NoQC_rmmissnp_${output}_g \
 --keep-allele-order \
 --make-bed \
 --allow-no-sex \
 --mind 0.2 \
 --out ${GWAS_path}/QC_${output}/QC/NoQC_rmmissnp_${output}_g_m

plink --bfile ${GWAS_path}/QC_${output}/QC/NoQC_rmmissnp_${output}_g_m \
 --keep-allele-order \
 --make-bed \
 --allow-no-sex \
 --maf 0.01 \
 --out ${GWAS_path}/QC_${output}/QC/NoQC_rmmissnp_${output}_g_m_maf 

plink --bfile ${GWAS_path}/QC_${output}/QC/NoQC_rmmissnp_${output}_g_m_maf \
 --keep-allele-order \
 --make-bed \
 --allow-no-sex \
 --hwe 1e-6 \
 --out ${GWAS_path}/QC_${output}/QC/QC_rmmissnp_${output}_g_m_maf_hwe

rm ${GWAS_path}/QC_${output}/MaMi/QC_rmmissnp_${Sample1}*
rm ${GWAS_path}/QC_${output}/MaMi/QC_rmmissnp_${Sample2}*
mv ${GWAS_path}/QC_${output}/QC/QC_rmmissnp_${output}_g_m_maf_hwe.bed ${GWAS_path}/QC_${output}/MaMi/rmmissnp_QC_${output}.bed
mv ${GWAS_path}/QC_${output}/QC/QC_rmmissnp_${output}_g_m_maf_hwe.bim ${GWAS_path}/QC_${output}/MaMi/rmmissnp_QC_${output}.bim
mv ${GWAS_path}/QC_${output}/QC/QC_rmmissnp_${output}_g_m_maf_hwe.fam ${GWAS_path}/QC_${output}/MaMi/rmmissnp_QC_${output}.fam
