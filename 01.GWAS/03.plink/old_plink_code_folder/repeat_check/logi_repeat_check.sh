#!/bin/sh
A=$1
sample_path="/scratch/x1997a11/GWAS/pdxen_AD/result_folder/BS_AD/QC_Imputed_AD_14_Imputed_AD_185/MaMi"
keep_path="/scratch/x1997a11/GWAS/pdxen_AD/Sample_folder/repeat_check"
result_path="/scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check"
inversion="/scratch/x1997a11/GWAS/pdxen_AD/reference_folder/inversion.txt"
mkdir ${result_path}/Case_Control/Case_Control${A}

plink --bfile ${sample_path}/QC_Imputed_AD_14_Imputed_AD_185_NoNA \
 --keep ${result_path}/Case/Case${A}/Case${A}.txt \
 --make-bed \
 --out ${result_path}/Case/Case${A}/raw_Case${A}

plink --bfile ${result_path}/Case/Case${A}/raw_Case${A} \
 --geno 0.2 \
 --make-bed \
 --out ${result_path}/Case/Case${A}/raw_Case${A}_g

plink --bfile ${result_path}/Case/Case${A}/raw_Case${A}_g \
 --mind 0.2 \
 --make-bed \
 --out ${result_path}/Case/Case${A}/raw_Case${A}_g_m

plink --bfile ${result_path}/Case/Case${A}/raw_Case${A}_g_m \
 --maf 0.01 \
 --make-bed \
 --out ${result_path}/Case/Case${A}/raw_Case${A}_g_m_maf

plink --bfile ${result_path}/Case/Case${A}/raw_Case${A}_g_m_maf \
 --hwe 1e-6 \
 --make-bed \
 --out ${result_path}/Case/Case${A}/Case${A}

plink --bfile ${sample_path}/QC_Imputed_AD_14_Imputed_AD_185_NoNA \
 --keep ${result_path}/Control/Control${A}/Control${A}.txt \
 --make-bed \
 --out ${result_path}/Control/Control${A}/raw_Control${A}

plink --bfile ${result_path}/Control/Control${A}/raw_Control${A} \
 --geno 0.2 \
 --make-bed \
 --out ${result_path}/Control/Control${A}/raw_Control${A}_g

plink --bfile ${result_path}/Control/Control${A}/raw_Control${A}_g \
 --mind 0.2 \
 --make-bed \
 --out ${result_path}/Control/Control${A}/raw_Control${A}_g_m

plink --bfile ${result_path}/Control/Control${A}/raw_Control${A}_g_m \
 --maf 0.01 \
 --make-bed \
 --out ${result_path}/Control/Control${A}/raw_Control${A}_g_m_maf

plink --bfile ${result_path}/Control/Control${A}/raw_Control${A}_g_m_maf \
 --hwe 1e-6 \
 --make-bed \
 --out ${result_path}/Control/Control${A}/Control${A}


#merge
plink --bfile ${result_path}/Case/Case${A}/Case${A} \
 --bmerge ${result_path}/Control/Control${A}/Control${A}.bed ${result_path}/Control/Control${A}/Control${A}.bim ${result_path}/Control/Control${A}/Control${A}.fam \
 --out ${result_path}/Case_Control/Case_Control${A}/NoQC_Case_Control${A}

#QC
plink --bfile ${result_path}/Case_Control/Case_Control${A}/NoQC_Case_Control${A} \
 --maf 0.01 \
 --make-bed \
 --out ${result_path}/Case_Control/Case_Control${A}/Case_Control${A}


#extract
plink --bfile ${result_path}/Case_Control/Case_Control${A}/Case_Control${A} \
 --assoc \
 --out ${result_path}/Case_Control/Case_Control${A}/Case_Control${A}_assoc

awk '$5 != "NA" && $6 != "NA" {print $0}' ${result_path}/Case_Control/Case_Control${A}/Case_Control${A}_assoc.assoc | grep -v "CHR" | awk '{print $2}' > ${result_path}/Case_Control/Case_Control${A}/extract_SNP_list.txt

plink --bfile ${result_path}/Case_Control/Case_Control${A}/Case_Control${A} \
 --extract ${result_path}/Case_Control/Case_Control${A}/extract_SNP_list.txt \
 --make-bed \
 --out ${result_path}/Case_Control/Case_Control${A}/raw_Case_Control${A}

#prune
plink --bfile ${result_path}/Case_Control/Case_Control${A}/raw_Case_Control${A} \
 --exclude ${inversion} \
 --range \
 --indep-pairwise 50 5 0.2 \
 --out ${result_path}/Case_Control/Case_Control${A}/raw_indepSNP

#genome
plink --bfile ${result_path}/Case_Control/Case_Control${A}/raw_Case_Control${A} \
 --extract ${result_path}/Case_Control/Case_Control${A}/raw_indepSNP.prune.in \
 --genome \
 --out ${result_path}/Case_Control/Case_Control${A}/raw_Case_Control${A}_genome

#MDS
mkdir ${result_path}/Case_Control/Case_Control${A}/MDS
plink --bfile ${result_path}/Case_Control/Case_Control${A}/raw_Case_Control${A} \
 --extract ${result_path}/Case_Control/Case_Control${A}/raw_indepSNP.prune.in \
 --make-bed \
 --out ${result_path}/Case_Control/Case_Control${A}/raw_indep_Case_Control${A}

plink --bfile ${result_path}/Case_Control/Case_Control${A}/raw_indep_Case_Control${A} \
 --read-genome ${result_path}/Case_Control/Case_Control${A}/raw_Case_Control${A}_genome.genome \
 --cluster --mds-plot 10 \
 --out ${result_path}/Case_Control/Case_Control${A}/MDS/Case_Control${A}_genome_MDS

mkdir ${result_path}/Case_Control/Case_Control${A}/PCA
plink --bfile ${result_path}/Case_Control/Case_Control${A}/raw_indep_Case_Control${A} \
 --double-id \
 --pca 10 \
 --set-missing-var-ids @:# \
 --out ${result_path}/Case_Control/Case_Control${A}/PCA/Case_Control${A}_PCA
awk '{print $1, $3, $4, $5, $6}' ${result_path}/Case_Control/Case_Control${A}/PCA/Case_Control${A}_PCA.eigenvec > ${result_path}/Case_Control/Case_Control${A}/PCA/raw_Case_Control${A}_PCA.csv
sed -i '1i\name PC1 PC2 PC3 PC4' ${result_path}/Case_Control/Case_Control${A}/PCA/raw_Case_Control${A}_PCA.csv

#logistic
mkdir ${result_path}/Case_Control/Case_Control${A}/logistic
#covar
awk '{print $1, $2, $4, $5, $6, $7 ,$8 ,$9 ,$10 ,$11, $12, $13}' ${result_path}/Case_Control/Case_Control${A}/MDS/Case_Control${A}_genome_MDS.mds > ${result_path}/Case_Control/Case_Control${A}/logistic.covar_mds.txt

plink --bfile ${result_path}/Case_Control/Case_Control${A}/raw_Case_Control${A} \
 --covar ${result_path}/Case_Control/Case_Control${A}/logistic.covar_mds.txt \
 --sex \
 --logistic \
 --hide-covar \
 --ci 0.95 \
 --out ${result_path}/Case_Control/Case_Control${A}/logistic/raw_Case_Control${A}_assoc

awk '!/'NA'/' ${result_path}/Case_Control/Case_Control${A}/logistic/raw_Case_Control${A}_assoc.assoc.logistic > ${result_path}/Case_Control/Case_Control${A}/logistic/Case_Control${A}_assoc.assoc.logistic






