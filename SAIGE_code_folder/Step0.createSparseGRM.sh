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
p2_Dir="/mnt/nas/gylee/Singularity/plink2"
#=================================================================
ana_dir="${Dir}/2.plink_result/NoImputed_${Sample}"
SAGE_CD="/mnt/nas/gylee/0.GWAS/Code/SAIGE_code_folder"

echo "=========================="
echo "Start Step0.createSparseGRM.sh"
echo "=========================="
mkdir ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LOG
Log="${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LOG/"

#Made by prune data
for A in $(seq 1 22)
do
nohup Rscript ${SAGE_CD}/createSparseGRM.R \
     --plinkFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC/chr${A}.LD_pruned_QC_${Sample} \
     --nThreads=32  \
     --outputPrefix=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${A}.sparseGRM \
     --numRandomMarkerforSparseKin=2000 \
     --relatednessCutoff=0.125 > ${Log}/Stop0.chr${A}.log 2>&1 &
done
