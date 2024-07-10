#!/bin/sh
if [ $# -ne 2 ];then
        echo "Please enter Sample_Name"
               exit
fi
Case=$1
Control=$2

GWAS_path="/scratch/x1997a11/GWAS/pdxen_AD/result_folder"
code_path="/scratch/x1997a11/GWAS/pdxen_AD/Code_folder"
Case_folder="${GWAS_path}/QC_${Case}"
Control_folder="${GWAS_path}/QC_${Control}"
#분석결과를 담을 파일
mkdir ${GWAS_path}/${Case}_${Control}_analysis_folder
analysis_folder="/${GWAS_path}/${Case}_${Control}_analysis_folder"
Case_pheno="${analysis_folder}/merge/Case_pheno.txt"
Control_pheno="${analysis_folder}/merge/Control_pheno.txt"
inversion="/scratch/x1997a11/GWAS/pdxen_AD/reference_folder/inversion.txt"
#=================================================================
echo "2_merge_CaseControl_analysis_major_missnp.sh" > ${analysis_folder}/merge/Used_2_merge_CaseControl_analysis_major_missnp.code
#================ missnp list ====================
echo "==================================="
echo "
		Missnp list 
"
echo "==================================="

mkdir ${analysis_folder}/merge

plink --bfile ${Case_folder}/MaMi/QC_${Case} \
 --bmerge ${Control_folder}/MaMi/QC_${Control}.bed ${Control_folder}/MaMi/QC_${Control}.bim ${Control_folder}/MaMi/QC_${Control}.fam \
 --out ${analysis_folder}/merge/raw_MaMi_${Case}_MaMi_${Control}

#missnp 이 없는 경우
awk '{print $2}' ${analysis_folder}/merge/raw_MaMi_${Case}_MaMi_${Control}.bim > ${analysis_folder}/merge/MaMi_${Case}_MaMi_${Control}.missnp
mv ${analysis_folder}/merge/raw_MaMi_${Case}_MaMi_${Control}.missnp ${analysis_folder}/merge/MaMi_${Case}_MaMi_${Control}.missnp
#================ missnp 제거=====================
echo "==================================="
echo "
                Missnp 제거
"
echo "==================================="

plink --bfile ${Case_folder}/MaMi/QC_${Case} \
 --exclude ${analysis_folder}/merge/MaMi_${Case}_MaMi_${Control}.missnp \
 --make-bed \
 --out ${analysis_folder}/merge/raw_NoQC_rmMissnp_${Case}

plink --bfile ${analysis_folder}/merge/raw_NoQC_rmMissnp_${Case} \
 --maf 0.01 \
 --make-bed \
 --out ${analysis_folder}/merge/rmMissnp_${Case}

plink --bfile ${Control_folder}/MaMi/QC_${Control} \
 --exclude ${analysis_folder}/merge/MaMi_${Case}_MaMi_${Control}.missnp \
 --make-bed \
 --out ${analysis_folder}/merge/raw_NoQC_rmMissnp_${Control}

plink --bfile ${analysis_folder}/merge/raw_NoQC_rmMissnp_${Control} \
 --maf 0.01 \
 --make-bed \
 --out ${analysis_folder}/merge/rmMissnp_${Control}


#=====================MaMi-MaMi merge=========================
plink --bfile  ${analysis_folder}/merge/rmMissnp_${Case} \
 --bmerge ${analysis_folder}/merge/rmMissnp_${Control}.bed ${analysis_folder}/merge/rmMissnp_${Control}.bim ${analysis_folder}/merge/rmMissnp_${Control}.fam \
 --out ${analysis_folder}/merge/raw_NoQC_MaMi_${Case}_MaMi_${Control}

echo "==================================="
echo "
                MAF QC ..ing
"
echo "==================================="
plink --bfile ${analysis_folder}/merge/raw_NoQC_MaMi_${Case}_MaMi_${Control} \
 --maf 0.01 \
 --make-bed \
 --out ${analysis_folder}/merge/raw_MaMi_${Case}_MaMi_${Control} 


echo "==================================="
echo "
                MAF QC ..end
"
echo "==================================="




plink --bfile ${analysis_folder}/merge/raw_MaMi_${Case}_MaMi_${Control} \
 --assoc \
 --out ${analysis_folder}/merge/raw_MaMi_${Case}_MaMi_${Control}_assoc

#grep -v "NA" ${analysis_folder}/merge/raw_MaMi_${Case}_MaMi_${Control}_assoc.assoc | grep -v "CHR" | awk '{print $2}' > ${analysis_folder}/merge/extract_SNP_list.txt
awk '$5 != "NA" && $6 != "NA" {print $0}' ${analysis_folder}/merge/raw_MaMi_${Case}_MaMi_${Control}_assoc.assoc | grep -v "CHR" | awk '{print $2}' > ${analysis_folder}/merge/extract_SNP_list.txt

#==========================================================

plink --bfile ${analysis_folder}/merge/raw_MaMi_${Case}_MaMi_${Control} \
 --extract ${analysis_folder}/merge/extract_SNP_list.txt \
 --make-bed \
 --out ${analysis_folder}/merge/raw_${Case}_${Control}_NoNA
#==========================================================

#prune
plink --bfile ${analysis_folder}/merge/raw_${Case}_${Control}_NoNA \
 --exclude ${inversion} \
 --range \
 --indep-pairwise 50 5 0.2 \
 --out ${analysis_folder}/merge/raw_indepSNP

#assoc
plink --bfile ${analysis_folder}/merge/raw_${Case}_${Control}_NoNA \
 --assoc \
 --out ${analysis_folder}/merge/raw_${Case}_${Control}_NoNA_assoc

awk '$5 != "NA" && $6 != "NA" {print $0}' ${analysis_folder}/merge/raw_${Case}_${Control}_NoNA_assoc.assoc > ${analysis_folder}/merge/${Case}_${Control}_NoNA_assoc.assoc

#genome
plink --bfile ${analysis_folder}/merge/raw_${Case}_${Control}_NoNA \
 --extract ${analysis_folder}/merge/raw_indepSNP.prune.in \
 --genome \
 --out ${analysis_folder}/merge/raw_${Case}_${Control}_NoNA_genome

#MDS
mkdir ${analysis_folder}/merge/MDS
plink --bfile ${analysis_folder}/merge/raw_${Case}_${Control}_NoNA \
 --extract ${analysis_folder}/merge/raw_indepSNP.prune.in \
 --make-bed \
 --out ${analysis_folder}/merge/raw_indep_${Case}_${Control}_NoNA

awk '{print $1, $2, $6}' ${analysis_folder}/merge/raw_${Case}_${Control}_NoNA.fam > ${analysis_folder}/merge/raw_pheno.txt
awk '$3 > 1 {print $0}' ${analysis_folder}/merge/raw_pheno.txt > ${analysis_folder}/merge/Case_pheno.txt
awk '$3 < 2 && $3 > 0 {print $0}' ${analysis_folder}/merge/raw_pheno.txt > ${analysis_folder}/merge/Control_pheno.txt 

plink --bfile ${analysis_folder}/merge/raw_indep_${Case}_${Control}_NoNA \
 --read-genome ${analysis_folder}/merge/raw_${Case}_${Control}_NoNA_genome.genome \
 --cluster --mds-plot 10 \
 --out ${analysis_folder}/merge/MDS/${Case}_${Control}_NoNA_genome_MDS

cat ${Case_pheno} ${Control_pheno} > ${analysis_folder}/merge/MDS/${Case}_${Control}_pheno.txt
awk '{print $1, $2, "Control"}' ${analysis_folder}/merge/MDS/${Case}_${Control}_pheno.txt > ${analysis_folder}/merge/MDS/Control_pheno.txt
awk '{print $1, $2, "Case"}' ${analysis_folder}/merge/MDS/${Case}_${Control}_pheno.txt > ${analysis_folder}/merge/MDS/Case_pheno.txt
cat ${analysis_folder}/merge/MDS/Control_pheno.txt ${analysis_folder}/merge/MDS/Case_pheno.txt | sed -e '1i\FID IID pheno' > ${analysis_folder}/merge/MDS/phenofile.txt
#Manhattan - QQplot code
cp /scratch/x1997a11/GWAS/pdxen_AD/Rcode/Manhattan_plot.R ${analysis_folder}/merge/
#PCA
mkdir ${analysis_folder}/merge/PCA
plink --bfile ${analysis_folder}/merge/raw_indep_${Case}_${Control}_NoNA \
 --double-id \
 --pca 10 \
 --set-missing-var-ids @:# \
 --out ${analysis_folder}/merge/PCA/${Case}_${Control}_NoNA_PCA

awk '{print $1, $3, $4, $5, $6}' ${analysis_folder}/merge/PCA/${Case}_${Control}_NoNA_PCA.eigenvec > ${analysis_folder}/merge/PCA/raw_${Case}_${Control}_PCA.csv
sed -i '1i\name PC1 PC2 PC3 PC4' ${analysis_folder}/merge/PCA/raw_${Case}_${Control}_PCA.csv

#logistic regression
mkdir ${analysis_folder}/merge/logistic
awk '{print $1, $2, $4, $5, $6, $7 ,$8 ,$9 ,$10 ,$11, $12, $13}' ${analysis_folder}/merge/MDS/${Case}_${Control}_NoNA_genome_MDS.mds > ${analysis_folder}/merge/logistic.covar_mds.txt
awk '{print $1, $5}' ${analysis_folder}/merge/raw_${Case}_${Control}_NoNA.fam > ${analysis_folder}/merge/logistic.covar_sex.txt
plink --bfile ${analysis_folder}/merge/raw_${Case}_${Control}_NoNA \
 --covar ${analysis_folder}/merge/logistic.covar_mds.txt \
 --sex \
 --logistic \
 --hide-covar \
 --ci 0.95 \
 --out ${analysis_folder}/merge/logistic/raw_${Case}_${Control}_NoNA_assoc

awk '!/'NA'/' ${analysis_folder}/merge/logistic/raw_${Case}_${Control}_NoNA_assoc.assoc.logistic > ${analysis_folder}/merge/logistic/${Case}_${Control}_NoNA_assoc.assoc.logistic

#haploview folder
mkdir ${analysis_folder}/merge/haploview

sed -i 's/    / /g' ${analysis_folder}/merge/logistic/${Case}_${Control}_NoNA_assoc.assoc.logistic
sed -i 's/    / /g' ${analysis_folder}/merge/${Case}_${Control}_NoNA_assoc.assoc
sed -i 's/    / /g' ${analysis_folder}/merge/logistic/${Case}_${Control}_NoNA_assoc.assoc.logistic
sed -i 's/    / /g' ${analysis_folder}/merge/${Case}_${Control}_NoNA_assoc.assoc
sed -i 's/  / /g' ${analysis_folder}/merge/logistic/${Case}_${Control}_NoNA_assoc.assoc.logistic
sed -i 's/  / /g' ${analysis_folder}/merge/${Case}_${Control}_NoNA_assoc.assoc
sed -i 's/  / /g' ${analysis_folder}/merge/logistic/${Case}_${Control}_NoNA_assoc.assoc.logistic
sed -i 's/  / /g' ${analysis_folder}/merge/${Case}_${Control}_NoNA_assoc.assoc
sed -i 's/  / /g' ${analysis_folder}/merge/logistic/${Case}_${Control}_NoNA_assoc.assoc.logistic
sed -i 's/  / /g' ${analysis_folder}/merge/${Case}_${Control}_NoNA_assoc.assoc

awk '$12 < 0.0005 {print $0}' ${analysis_folder}/merge/logistic/${Case}_${Control}_NoNA_assoc.assoc.logistic > ${analysis_folder}/merge/logistic/low_${Case}_${Control}_NoNA_assoc.assoc.logistic
sed -i '1i\ CHR SNP BP A1 TEST NMISS OR SE L95 U95 STAT P' ${analysis_folder}/merge/logistic/low_${Case}_${Control}_NoNA_assoc.assoc.logistic
awk '$9 < 0.0005 {print $0}' ${analysis_folder}/merge/${Case}_${Control}_NoNA_assoc.assoc > ${analysis_folder}/merge/low_${Case}_${Control}_NoNA_assoc.assoc
sed -i '1i\ CHR SNP BP A1 F_A F_U A2 CHISQ P OR' ${analysis_folder}/merge/low_${Case}_${Control}_NoNA_assoc.assoc

mv ${analysis_folder}/merge/raw_${Case}_${Control}_NoNA.bed ${analysis_folder}/merge/${Case}_${Control}_NoNA.bed
mv ${analysis_folder}/merge/raw_${Case}_${Control}_NoNA.bim ${analysis_folder}/merge/${Case}_${Control}_NoNA.bim
mv ${analysis_folder}/merge/raw_${Case}_${Control}_NoNA.fam ${analysis_folder}/merge/${Case}_${Control}_NoNA.fam
mv ${analysis_folder}/merge/raw_${Case}_${Control}_NoNA.log ${analysis_folder}/merge/${Case}_${Control}_NoNA.log
rm ${analysis_folder}/merge/raw*

