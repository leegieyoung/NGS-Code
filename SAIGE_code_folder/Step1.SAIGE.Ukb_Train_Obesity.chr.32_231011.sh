#!/bin/sh
if [ $# -ne 3 ];then
        echo "Please enter Sample_Name, thread(s)"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
#beagle 진행시Same position 제거를 안했을 경우 사용하는 코드
Sample=$1
Thread=$2
Chr=$3
Dir="/mnt/nas/gylee/0.GWAS"
impute_dir="${Dir}/1.Input/${Sample}/0.Impute"
QC_dir="${Dir}/1.Input/${Sample}/1.QC"
MaMi_dir="${Dir}/1.Input/${Sample}/2.MaMi"
Case_pheno="${ana_dir}/temp/11.Case_pheno.txt"
Control_pheno="${ana_dir}/temp/11.Control_pheno.txt"
inversion="${Dir}/REFERENCE/inversion.txt"
Rcode="/mnt/nas/gylee/0.GWAS/Code/R"
p2_Dir="/mnt/nas/gylee/Singurality/plink2"
Output_dir="${Dir}/2.plink_result/${Sample}"
#=================================================================
ana_dir="${Dir}/2.plink_result/NoImputed_${Sample}"
SAIGE_CD="/mnt/nas/gylee/0.GWAS/Code/SAIGE_code_folder"

#mkdir ${Dir}/2.plink_result/${Sample}/chr${Chr}
#cp ${Dir}/2.plink_result/${Sample}/chr${Chr}.imputation_g_m_maf_hwe_bfile.bim ${Dir}/2.plink_result/${Sample}/chr${Chr}/

echo "=========================="
echo "Start Step0.createSparseGRM.sh"
echo "=========================="
mkdir -p ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LOG
mkdir -p ${Dir}/2.SAIGE_result/${Sample}/result

Log="${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LOG"

export TRAIT1="${Sample}"
export CHR="${Chr}"
Rscript ${SAIGE_CD}/SAIGE.Ukb_Train_Obesity.chr21_22-32.R

for POS in $(seq 1 32)
do
/mnt/nas/gylee/Singurality/plink2/plink2 --pfile ${Dir}/2.plink_result/${Sample}/chr${Chr}.imputation_g_m_maf_hwe \
 --extract ${Dir}/2.plink_result/${Sample}/chr${Chr}/chr${Chr}-${POS}.list \
 --make-bed \
 --out ${Dir}/2.plink_result/${Sample}/chr${Chr}/chr${Chr}-${POS}
done
