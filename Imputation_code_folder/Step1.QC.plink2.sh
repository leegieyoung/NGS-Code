#!/bin/sh
if [ $# -ne 1 ];then
        echo "Please enter Input file name"
               exit
fi
Input=$1
Dir="/ichrogene/project/temp/gylee/0.GWAS/05.Imputation" #Change to your path
Image="/ichrogene/project/temp/gylee/Singularity"  #Change to your image file path
p2_Dir="/ichrogene/project/temp/gylee/Singularity/plink2"
p19_Dir="/ichrogene/project/temp/gylee/Singularity/plink19"
mkdir ${Dir}/temp
mkdir ${Dir}/1.shapeit #Result folder of shapit4 that phasing tool
mkdir ${Dir}/2.minimac4  #Result folder of minimac4 that impuation tool

#ped to 
# 1. Remove N and Make list of duplicated variant

awk '{if($5!="N") print $2}' ${Dir}/0.INPUT/${Input}/${Input}.bim > ${Dir}/0.INPUT/${Input}/${Input}.rmN

${p19_Dir}/plink --bfile ${Dir}/0.INPUT/${Input}/${Input} \
 --extract ${Dir}/0.INPUT/${Input}/${Input}.rmN \
 --keep-allele-order \
 --allow-no-sex \
 --make-bed \
 --out ${Dir}/0.INPUT/${Input}/${Input}.rmN

${p2_Dir}/plink2 --bfile ${Dir}/0.INPUT/${Input}/${Input}.rmN \
 --make-pgen \
 --out ${Dir}/0.INPUT/${Input}/${Input}
 
${p2_Dir}/plink2 --pfile ${Dir}/0.INPUT/${Input}/${Input} \
 --make-just-bim \
 --rm-dup force-first \
 --out ${Dir}/temp/rmdupli


# 2. Remove INDEL
awk '$5 != "-" && $6 != "-" && $5 != "I" && $6 != "I" && $5 != "D" && $6 != "D" && $5 != "0" && $6 != "0" && $5 != "N" && $6 != "N" && $5 != "." && $6 != "." ' ${Dir}/temp/rmdupli.bim > ${Dir}/temp/rmdupli_rmINDEL.bim

awk '{print $2}' ${Dir}/temp/rmdupli_rmINDEL.bim > ${Dir}/temp/rmdupli_rmINDEL.list

# 3. Only acgt & maf 0.01
${p2_Dir}/plink2 --pfile ${Dir}/0.INPUT/${Input}/${Input} \
 --snps-only just-acgt \
 --extract ${Dir}/temp/rmdupli_rmINDEL.list \
 --recode vcf \
 --make-pgen \
 --out ${Dir}/temp/onlysnp

singularity exec --bind /ichrogene/:/ichrogene/ /ichrogene/project/temp/gylee/Singularity/minimac4_1.0.3.sif bcftools norm -d all -O z --threads 30 ${Dir}/temp/onlysnp.vcf -o ${Dir}/temp/onlysnp.norm.vcf

grep -v '^#'  ${Dir}/temp/onlysnp.norm.vcf | awk '{print $3}' > ${Dir}/temp/onlysnp.norm.list

${p2_Dir}/plink2 --pfile ${Dir}/temp/onlysnp \
 --extract ${Dir}/temp/onlysnp.norm.list \
 --geno 0.2 \
 --make-pgen \
 --out ${Dir}/temp/onlysnp_g

${p2_Dir}/plink2 --pfile ${Dir}/temp/onlysnp_g \
 --mind 0.2 \
 --make-pgen \
 --out ${Dir}/temp/onlysnp_g_m

${p2_Dir}/plink2 --pfile ${Dir}/temp/onlysnp_g_m \
 --maf 0.01 \
 --make-pgen \
 --out ${Dir}/0.INPUT/${Input}/imputation_${Input}_rmdupli_rmINDEL

rm -rf ${Dir}/temp

for Chr in $(seq 1 22)
do
${p2_Dir}/plink2 --pfile ${Dir}/0.INPUT/${Input}/imputation_${Input}_rmdupli_rmINDEL \
        --chr ${Chr} \
        --recode vcf bgz \
        --out ${Dir}/0.INPUT/${Input}/imputation_${Input}_${Chr}_rmdupli_rmINDEL
done
