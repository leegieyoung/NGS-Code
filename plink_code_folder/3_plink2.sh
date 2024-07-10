#!/bin/sh
if [ $# -ne 2 ];then
        echo "Please enter Sample_Name, thread(s)"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
#beagle 진행시Same position 제거를 안했을 경우 사용하는 코드
Sample=$1
Thread=$2
Dir="/mnt/nas/gylee/0.GWAS"
impute_dir="${Dir}/1.Input/${Sample}/0.Impute"
QC_dir="${Dir}/1.Input/${Sample}/1.QC"
MaMi_dir="${Dir}/1.Input/${Sample}/2.MaMi"
Case_pheno="${ana_dir}/temp/11.Case_pheno.txt"
Control_pheno="${ana_dir}/temp/11.Control_pheno.txt"
inversion="${Dir}/REFERENCE/inversion.txt"
Rcode="/mnt/nas/gylee/0.GWAS/Code/R"
p2_Dir="/mnt/nas/gylee/Singurality/plink2"
#=================================================================

mkdir -p ${Dir}/2.plink_result/NoImputed_${Sample}/temp
mkdir -p ${Dir}/2.plink_result/NoImputed_${Sample}/result
ana_dir="${Dir}/2.plink_result/NoImputed_${Sample}"

#=========What is Code ? =======================================
echo "${code_path}/2_Single_version" > ${ana_dir}/2_single_version

#PCA
mkdir ${ana_dir}/PCA
${p2_Dir}/plink2 --bfile ${ana_dir}/temp/8.raw_${Sample}_NoNA \
 --double-id \
 --pca 20 approx \
 --memory 240000 \
 --set-missing-var-ids @:# \
 --threads ${Thread} \
 --out ${ana_dir}/PCA/${Sample}_NoNA_PCA

awk '{print $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22}' ${ana_dir}/PCA/${Sample}_NoNA_PCA.eigenvec > ${ana_dir}/PCA/raw_${Sample}_PCA.csv

sed -i '1i\name PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 PC11 PC12 PC13 PC14 PC15 PC16 PC17 PC18 PC19 PC20' ${ana_dir}/PCA/raw_${Sample}_PCA.csv
paste -d '\t' ${ana_dir}/PCA/raw_${Sample}_PCA.csv ${ana_dir}/PCA/${Sample}_NoNA_PCA.eigenval > ${ana_dir}/PCA/raw_${Sample}_PCA-eigenval.csv

awk '{print $6}' ${ana_dir}/temp/8.raw_${Sample}_NoNA.fam > ${ana_dir}/PCA/${Sample}_NoNA.pheno
sed -i -e 's/2/Case/g' -e 's/1/Control/g' ${ana_dir}/PCA/${Sample}_NoNA.pheno
sed -i '1iname' ${ana_dir}/PCA/${Sample}_NoNA.pheno
awk '{print $1}' ${ana_dir}/temp/8.raw_${Sample}_NoNA.fam > ${ana_dir}/PCA/${Sample}_NoNA.sample
sed -i '1isample' ${ana_dir}/PCA/${Sample}_NoNA.sample
awk '{print $3, $4}' ${ana_dir}/PCA/${Sample}_NoNA_PCA.eigenvec > ${ana_dir}/PCA/raw_${Sample}_PCA.PC12
#sed -i '1iPC1 PC2' ${ana_dir}/PCA/raw_${Sample}_PCA.PC12
paste -d ' ' ${ana_dir}/PCA/${Sample}_NoNA.pheno ${ana_dir}/PCA/raw_${Sample}_PCA.PC12 ${ana_dir}/PCA/${Sample}_NoNA.sample > ${ana_dir}/PCA/PCA.txt

## Add phenotype
export ANA_DIR="${Dir}/2.plink_result/NoImputed_${Sample}/"
export DIR="${Dir}/2.SAIGE_result/${Sample}/sparseGRM/"
export SAMPLE="${Sample}"

mkdir -p ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/
Rscript /mnt/nas/gylee/0.GWAS/Code/plink_code_folder/4_test.R

awk '{print $1,$1}' ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/${Sample}.covariate.txt > ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/${Sample}.covariate.IID

${p2_Dir}/plink2 --bfile ${ana_dir}/temp/8.raw_${Sample}_NoNA \
 --keep ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/${Sample}.covariate.IID \
 --make-bed \
 --out ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC_${Sample}

for A in $(seq 1 22)
do
${p2_Dir}/plink2 --bfile ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC_${Sample} \
--chr ${A} \
 --make-bed \
 --out ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${A}.LD_pruned_QC_${Sample}
done

