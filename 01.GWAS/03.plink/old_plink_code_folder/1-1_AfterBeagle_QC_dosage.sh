#!/bin/sh
if [ $# -ne 1 ];then
        echo "Please enter Sample_Name"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
#beagle 진행시Same position 제거를 안했을 경우 사용하는 코드
Sample=$1
GWAS_path="/scratch/hpc46a05/GWAS/result"
code_path="/scratch/hpc46a05/GWAS/Code/plink_code_folder"
mkdir ${GWAS_path}/QC_Imputed_${Sample}
result_folder="${GWAS_path}/QC_Imputed_${Sample}"
#imputed_sample folder
beagle_folder="/scratch/hpc46a05/GWAS/beagle_result"
Sample_folder="${beagle_folder}/3.${Sample}_imputed_${Sample}"
raw_folder="/scratch/hpc46a05/GWAS/raw_data/"
#beagle_folder="/scratch/x1997a11/GWAS/pdxen_AD/beagle/Result_folder"

#1.SampleQC(geno mind impute-sex hwe)
mkdir ${result_folder}/raw_merge
mkdir ${result_folder}/QC
mkdir ${result_folder}/MaMi

cp -v ${beagle_folder}/imputed_QC_${Sample}/* ${Sample_folder}/
for Chrom in {1..22}
do
plink --vcf ${Sample_folder}/CHR${Chrom}_${Sample}_noIndel_rmDupID_beagle.vcf.gz \
 --double-id --make-bed --keep-allele-order \
 --out ${result_folder}/raw_merge/NoQC_CHR${Chrom}_${Sample}
done

rmIDdot
for Chrom in {1..22}
do
echo ''${Chrom}_'rmIDdot'
awk '$2 != "." && $5 != "." && $6 != "." {print $2}' ${result_folder}/raw_merge/NoQC_CHR${Chrom}_${Sample}.bim > ${result_folder}/raw_merge/NoQC_rmIDdot_CHR${Chrom}_${Sample}.txt
plink --bfile ${result_folder}/raw_merge/NoQC_CHR${Chrom}_${Sample} \
 --extract ${result_folder}/raw_merge/NoQC_rmIDdot_CHR${Chrom}_${Sample}.txt \
 --make-bed --keep-allele-order \
 --out ${result_folder}/raw_merge/NoQC_rmIDdot_CHR${Chrom}_${Sample}
done

for Chrom in {2..22}
do
echo ${result_folder}/raw_merge/NoQC_rmIDdot_CHR${Chrom}_${Sample}.bed ${result_folder}/raw_merge/NoQC_rmIDdot_CHR${Chrom}_${Sample}.bim ${result_folder}/raw_merge/NoQC_rmIDdot_CHR${Chrom}_${Sample}.fam > ${result_folder}/raw_merge/mergelist_${Chrom}.txt
done

cat ${result_folder}/raw_merge/mergelist_*.txt > ${result_folder}/raw_merge/mergelist.txt

plink --bfile ${result_folder}/raw_merge/NoQC_rmIDdot_CHR1_${Sample} \
 --merge-list ${result_folder}/raw_merge/mergelist.txt \
 --make-bed \
 --keep-allele-order \
 --out ${result_folder}/QC/raw_NoQC_same_posi_Imputed_${Sample} 

#rm Same position
grep '^Warning' ${result_folder}/QC/raw_NoQC_same_posi_Imputed_${Sample}.log > ${result_folder}/QC/same_position.log
sed -i "s/'/\t/g" ${result_folder}/QC/same_position.log
awk '{print $3}' ${result_folder}/QC/same_position.log > ${result_folder}/QC/same_position.rsid
plink --bfile ${result_folder}/QC/raw_NoQC_same_posi_Imputed_${Sample} \
 --exclude ${result_folder}/QC/same_position.rsid \
 --keep-allele-order \
 --make-bed \
 --out ${result_folder}/QC/raw_NoQC_Imputed_${Sample}

plink --bfile ${result_folder}/QC/raw_NoQC_Imputed_${Sample} \
 --geno 0.05 \
 --keep-allele-order \
 --make-bed \
 --out ${result_folder}/QC/raw_NoQC_Imputed_${Sample}_g

plink --bfile ${result_folder}/QC/raw_NoQC_Imputed_${Sample}_g \
 --mind 0.2 \
 --keep-allele-order \
 --make-bed \
 --out ${result_folder}/QC/raw_NoQC_Imputed_${Sample}_g_m

plink --bfile ${result_folder}/QC/raw_NoQC_Imputed_${Sample}_g_m \
 --maf 0.05 \
 --keep-allele-order \
 --make-bed \
 --out ${result_folder}/QC/raw_NoQC_Imputed_${Sample}_g_m_maf

plink --bfile ${result_folder}/QC/raw_NoQC_Imputed_${Sample}_g_m_maf \
 --hwe 1e-6 \
 --keep-allele-order \
 --make-bed \
 --out ${result_folder}/MaMi/QC_Imputed_${Sample}

mkdir ${result_folder}/dosage
echo '================================'
echo '                                '
echo '      Do making dosage file  '
echo '        Frq is all dataset '
echo '                                '
echo '================================'
plink --bfile ${result_folder}/MaMi/QC_Imputed_${Sample} \
 --keep-allele-order \
 --freq \
 --out ${result_folder}/dosage/raw_QC_Imputed_${Sample}_freq

awk '{print $1, $2, $3 ,$4, $5}' ${result_folder}/dosage/raw_QC_Imputed_${Sample}_freq.frq > ${result_folder}/dosage/QC_Imputed_${Sample}_freq.frq

plink --bfile ${result_folder}/MaMi/QC_Imputed_${Sample} \
 --recode vcf-iid \
 --keep-allele-order \
 --out ${result_folder}/dosage/raw_QC_Imputed_${Sample}

grep -v "^##" ${result_folder}/dosage/raw_QC_Imputed_${Sample}.vcf > ${result_folder}/dosage/QC_Imputed_${Sample}.vcf
cat ${result_folder}/dosage/QC_Imputed_${Sample}.vcf | colrm 1 9 > ${result_folder}/dosage/raw_QC_Imputed_${Sample}.dosage
sed -i -e 's/0\/0/0/g' -e 's/0\/1/1/g' -e 's/1\/1/2/g' ${result_folder}/dosage/raw_QC_Imputed_${Sample}.dosage
paste -d ' ' ${result_folder}/dosage/QC_Imputed_${Sample}_freq.frq ${result_folder}/dosage/raw_QC_Imputed_${Sample}.dosage > ${result_folder}/dosage/QC_Imputed_${Sample}.dosage

sed -n '1p' ${result_folder}/dosage/QC_Imputed_${Sample}.dosage > ${result_folder}/dosage/dosage.head
for A in $(seq 1 22)
do
grep -w "^${A}" ${result_folder}/dosage/QC_Imputed_${Sample}.dosage > ${result_folder}/dosage/raw_QC_Imputed_chr${A}_${Sample}.dosage
cat ${result_folder}/dosage/dosage.head ${result_folder}/dosage/raw_QC_Imputed_chr${A}_${Sample}.dosage > ${result_folder}/dosage/QC_Imputed_chr${A}_${Sample}.dosage
yes n | gzip ${result_folder}/dosage/QC_Imputed_chr${A}_${Sample}.dosage
done



echo ""
echo "===================Delete Used Data========================"
echo ""
rm ${result_folder}/raw_merge/*.bed
rm ${result_folder}/raw_merge/*.bim
rm ${result_folder}/raw_merge/*.fam
rm ${result_folder}/QC/*.bed
rm ${result_folder}/QC/*.bim
rm ${result_folder}/QC/*.fam
rm ${result_folder}/dosage/*.dosage
rm ${result_folder}/dosage/raw*
