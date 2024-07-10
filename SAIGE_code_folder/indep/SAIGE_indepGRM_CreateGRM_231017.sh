#!/bin/sh
if [ $# -ne 2 ];then
        echo "Please enter Sample_Name, thread(s)"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
#beagle 진행시Same position 제거를 안했을 경우 사용하는 코드
Sample=$1
Thread=$2
Dir="/ichrogene/project/temp/gylee/0.GWAS"
impute_dir="${Dir}/1.Input/${Sample}/0.Impute"
QC_dir="${Dir}/1.Input/${Sample}/1.QC"
inversion="${Dir}/REFERENCE/inversion.txt"
Rcode="/ichrogene/project/temp/gylee/0.GWAS/Code/R"
p2_Dir="/ichrogene/project/temp/gylee/Singularity/plink2"
#=================================================================
ana_dir="${Dir}/2.plink_result/NoImputed_${Sample}"
SAGE_CD="/ichrogene/project/temp/gylee/0.GWAS/Code/SAIGE_code_folder"

echo "=========================="
echo "Start Step0.createSparseGRM.sh"
echo "=========================="
mkdir ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LOG
Log="${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LOG/"
mkdir ${Log}/createSparseGRM
#Made by prune data
for A in $(seq 1 20)
do
Rscript ${SAGE_CD}/createSparseGRM.R \
     --plinkFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC/chr${A}.LD_pruned_QC_${Sample} \
     --nThreads=27  \
     --outputPrefix=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${A}.sparseGRM \
     --numRandomMarkerforSparseKin=2000 \
     --relatednessCutoff=0.125 > ${Log}/createSparseGRM/Stop0.chr${A}.log 
done

for A in $(seq 21 22)
do
Rscript ${SAGE_CD}/createSparseGRM.R \
     --plinkFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC/chr${A}.LD_pruned_QC_${Sample} \
     --nThreads=27  \
     --outputPrefix=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${A}.sparseGRM \
     --numRandomMarkerforSparseKin=1000 \
     --relatednessCutoff=0.125 > ${Log}/createSparseGRM/Stop0.chr${A}.log 
done

