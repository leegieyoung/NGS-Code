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
#=================================================================

mkdir -p ${Dir}/2.plink_result/NoImputed_${Sample}/temp
mkdir -p ${Dir}/2.plink_result/NoImputed_${Sample}/result
ana_dir="${Dir}/2.plink_result/NoImputed_${Sample}"

#=========What is Code ? =======================================
echo "${code_path}/2_Single_version" > ${ana_dir}/2_single_version

awk '{print $1, $2, $6}' ${ana_dir}/temp/8.raw_${Sample}_NoNA.fam > ${ana_dir}/temp/11.raw_pheno.txt
awk '$3 > 1 {print $0}' ${ana_dir}/temp/11.raw_pheno.txt > ${ana_dir}/temp/11.Case_pheno.txt
awk '$3 < 2 && $3 > 0 {print $0}' ${ana_dir}/temp/11.raw_pheno.txt > ${ana_dir}/temp/11.Control_pheno.txt 

plink --bfile ${ana_dir}/temp/11.raw_indep_${Sample}_NoNA \
 --read-genome ${ana_dir}/temp/10.raw_${Sample}_NoNA_genome.genome \
 --cluster --mds-plot 10 \
 --threads ${Thread} \
 --allow-no-sex \
 --out ${ana_dir}/MDS/${Sample}_NoNA_genome_MDS

cat ${Case_pheno} ${Control_pheno} > ${ana_dir}/MDS/${Sample}_pheno.txt
awk '{print $1, $2, "Control"}' ${ana_dir}/temp/11.Control_pheno.txt > ${ana_dir}/MDS/Control_pheno.txt
awk '{print $1, $2, "Case"}' ${ana_dir}/temp/11.Case_pheno.txt > ${ana_dir}/MDS/Case_pheno.txt
cat ${ana_dir}/MDS/Control_pheno.txt ${ana_dir}/MDS/Case_pheno.txt | sed -e '1i\FID IID pheno' > ${ana_dir}/MDS/phenofile.txt
cp ${ana_dir}/MDS/${Sample}_NoNA_genome_MDS.mds > ${ana_dir}/MDS/MDS.mds
cp ${Rcode}/MDS_plot.R ${ana_dir}/MDS/MDS_plot.R

#PCA
mkdir ${ana_dir}/PCA
plink --bfile ${ana_dir}/temp/8.raw_${Sample}_NoNA \
 --double-id \
 --pca 10 \
 --allow-no-sex \
 --set-missing-var-ids @:# \
 --threads ${Thread} \
 --out ${ana_dir}/PCA/${Sample}_NoNA_PCA

awk '{print $1, $3, $4, $5, $6}' ${ana_dir}/PCA/${Sample}_NoNA_PCA.eigenvec > ${ana_dir}/PCA/raw_${Sample}_PCA.csv
sed -i '1i\name PC1 PC2 PC3 PC4' ${ana_dir}/PCA/raw_${Sample}_PCA.csv
paste -d '\t' ${ana_dir}/PCA/raw_${Sample}_PCA.csv ${ana_dir}/PCA/${Sample}_NoNA_PCA.eigenval > ${ana_dir}/PCA/raw_${Sample}_PCA-eigenval.csv

awk '{print $6}' ${ana_dir}/temp/8.raw_${Sample}_NoNA.fam > ${ana_dir}/PCA/${Sample}_NoNA.pheno
sed -i -e 's/2/Case/g' -e 's/1/Control/g' ${ana_dir}/PCA/${Sample}_NoNA.pheno
sed -i '1iname' ${ana_dir}/PCA/${Sample}_NoNA.pheno
awk '{print $1}' ${ana_dir}/temp/8.raw_${Sample}_NoNA.fam > ${ana_dir}/PCA/${Sample}_NoNA.sample
sed -i '1isample' ${ana_dir}/PCA/${Sample}_NoNA.sample
awk '{print $3, $4}' ${ana_dir}/PCA/${Sample}_NoNA_PCA.eigenvec > ${ana_dir}/PCA/raw_${Sample}_PCA.PC12
sed -i '1iPC1 PC2' ${ana_dir}/PCA/raw_${Sample}_PCA.PC12
paste -d ' ' ${ana_dir}/PCA/${Sample}_NoNA.pheno ${ana_dir}/PCA/raw_${Sample}_PCA.PC12 ${ana_dir}/PCA/${Sample}_NoNA.sample > ${ana_dir}/PCA/PCA.txt

#logistic regression (Odd Ratio)
mkdir ${ana_dir}/logistic
awk '{print $1, $2, $4, $5, $6, $7 ,$8 ,$9 ,$10 ,$11, $12, $13}' ${ana_dir}/MDS/${Sample}_NoNA_genome_MDS.mds > ${ana_dir}/logistic.covar_mds.txt

plink --bfile ${ana_dir}/temp/8.raw_${Sample}_NoNA \
 --covar ${ana_dir}/logistic.covar_mds.txt \
 --logistic \
 --sex \
 --hide-covar \
 --ci 0.95 \
 --threads ${Thread} \
 --allow-no-sex \
 --out ${ana_dir}/logistic/raw_${Sample}_NoNA_assoc

#awk '!/'NA'/' ${ana_dir}/logistic/raw_${Sample}_NoNA_assoc.assoc.logistic > ${ana_dir}/logistic/${Sample}_NoNA_assoc.assoc.logistic
awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' ${ana_dir}/logistic/raw_${Sample}_NoNA_assoc.assoc.logistic > ${ana_dir}/logistic/${Sample}_NoNA_assoc.assoc.logistic

#Manhattan - QQplot code
cp ${Rcode}/Manhattan_plot.R ${ana_dir}/logistic/

#haploview folder
mkdir ${ana_dir}/haploview

mv ${ana_dir}/temp/8.raw_${Sample}_NoNA.bed ${ana_dir}/${Sample}_NoNA.bed
mv ${ana_dir}/temp/8.raw_${Sample}_NoNA.bim ${ana_dir}/${Sample}_NoNA.bim
mv ${ana_dir}/temp/8.raw_${Sample}_NoNA.fam ${ana_dir}/${Sample}_NoNA.fam
mv ${ana_dir}/temp/8.raw_${Sample}_NoNA.log ${ana_dir}/${Sample}_NoNA.log

#rm ${ana_dir}/raw*
mkdir ${ana_dir}/raw_file
mv ${ana_dir}/raw* ${ana_dir}/raw_file

#anno
mkdir ${ana_dir}/logistic/anno

plink --bfile ${ana_dir}/${Sample}_NoNA \
 --freq case-control \
 --threads ${Thread} \
 --allow-no-sex \
 --out ${ana_dir}/logistic/anno/raw_${Sample}_NoNA_freq

plink --bfile ${ana_dir}/${Sample}_NoNA \
 --freq \
 --threads ${Thread} \
 --allow-no-sex \
 --out ${ana_dir}/logistic/anno/${Sample}_NoNA_freq
 
awk '{print $1, $2, $3, $4, $5, $6, $7, $8}' ${ana_dir}/logistic/anno/raw_${Sample}_NoNA_freq.frq.cc > ${ana_dir}/logistic/anno/${Sample}_NoNA_freq.frq.cc

#anno_summary_file
#chr
awk '{print $1}'  ${ana_dir}/logistic/anno/${Sample}_NoNA_freq.frq.cc > ${ana_dir}/logistic/anno/${Sample}_NoNA.chr
#posi
awk '{print $3}' ${ana_dir}/logistic/raw_${Sample}_NoNA_assoc.assoc.logistic > ${ana_dir}/logistic/anno/${Sample}_NoNA.posi
#rsID
awk '{print $2}' ${ana_dir}/logistic/raw_${Sample}_NoNA_assoc.assoc.logistic > ${ana_dir}/logistic/anno/${Sample}_NoNA.rsID
#A1-A2
awk '{print $3, $4}' ${ana_dir}/logistic/anno/${Sample}_NoNA_freq.frq.cc > ${ana_dir}/logistic/anno/${Sample}_NoNA.A1-A2
awk '{print $3}' ${ana_dir}/logistic/anno/${Sample}_NoNA_freq.frq.cc > ${ana_dir}/logistic/anno/${Sample}_NoNA.A1
awk '{print $4}' ${ana_dir}/logistic/anno/${Sample}_NoNA_freq.frq.cc > ${ana_dir}/logistic/anno/${Sample}_NoNA.A2
#OR
awk '{print $7}' ${ana_dir}/logistic/raw_${Sample}_NoNA_assoc.assoc.logistic > ${ana_dir}/logistic/anno/${Sample}_NoNA.OR
#pval
awk '{print $12}' ${ana_dir}/logistic/raw_${Sample}_NoNA_assoc.assoc.logistic > ${ana_dir}/logistic/anno/${Sample}_NoNA.pval
#SE
awk '{print $8}' ${ana_dir}/logistic/raw_${Sample}_NoNA_assoc.assoc.logistic > ${ana_dir}/logistic/anno/${Sample}_NoNA.SE
#maf
awk '{print $5}' ${ana_dir}/logistic/anno/${Sample}_NoNA_freq.frq > ${ana_dir}/logistic/anno/${Sample}_NoNA.maf


echo '================================'
echo '                                '
echo '      Do making summary file    '
echo '                                '
echo '================================'

paste ${ana_dir}/logistic/anno/${Sample}_NoNA_freq.frq.cc \
 ${ana_dir}/logistic/anno/${Sample}_NoNA.posi \
 ${ana_dir}/logistic/anno/${Sample}_NoNA.pval \
 > ${ana_dir}/logistic/anno/summary_result.csv

echo '================================'
echo '                                '
echo ' Do making predixcan input file '
echo '                                '
echo '================================'
paste -d '\t' ${ana_dir}/logistic/anno/${Sample}_NoNA.rsID \
 ${ana_dir}/logistic/anno/${Sample}_NoNA.chr \
 ${ana_dir}/logistic/anno/${Sample}_NoNA.posi \
 ${ana_dir}/logistic/anno/${Sample}_NoNA.A1-A2 \
 ${ana_dir}/logistic/anno/${Sample}_NoNA.maf \
 ${ana_dir}/logistic/anno/${Sample}_NoNA.OR \
 ${ana_dir}/logistic/anno/${Sample}_NoNA.pval \
 > ${ana_dir}/logistic/anno/raw_${Sample}.txt

grep -v 'NA' ${ana_dir}/logistic/anno/raw_${Sample}.txt > ${ana_dir}/logistic/anno/${Sample}.txt
sed -i 's/ /\t/g' ${ana_dir}/logistic/anno/${Sample}.txt

yes n | gzip ${ana_dir}/logistic/anno/${Sample}.txt 

echo '================================'
echo '                                '
echo '    Do making FUMA input file '
echo '                                '
echo '================================'

mkdir ${ana_dir}/logistic/FUMA
paste -d '\t' ${ana_dir}/logistic/anno/${Sample}_NoNA.chr \
 ${ana_dir}/logistic/anno/${Sample}_NoNA.posi \
 ${ana_dir}/logistic/anno/${Sample}_NoNA.rsID \
 ${ana_dir}/logistic/anno/${Sample}_NoNA.pval \
 ${ana_dir}/logistic/anno/${Sample}_NoNA.A1 \
 ${ana_dir}/logistic/anno/${Sample}_NoNA.A2 \
 ${ana_dir}/logistic/anno/${Sample}_NoNA.OR \
 ${ana_dir}/logistic/anno/${Sample}_NoNA.SE > ${ana_dir}/logistic/FUMA/fuma_input.txt

sed -i '1d' ${ana_dir}/logistic/FUMA/fuma_input.txt
sed -i '1ichromosome\tposition\tSNP\tP-value\tA1\tA2\tOR\tSE' ${ana_dir}/logistic/FUMA/fuma_input.txt

echo ""
echo "===========================Delete Used Data==========================="
echo ""
rm ${ana_dir}/logistic/anno/*_freq.frq
rm ${ana_dir}/logistic/anno/*_freq.frq.cc
rm ${ana_dir}/logistic/anno/${Sample}.txt
rm ${ana_dir}/logistic/raw_${Sample}_NoNA_assoc.assoc.logistic
#rm ${ana_dir}/raw_file/raw_${Sample}_assoc.assoc

