#!/bin/sh
if [ $# -ne 2 ];then
        echo "Please enter Sample_Name, thread(s)"
               exit
fi
Sample=$1
Thread=$2
Dir="/ichrogene/project/temp/gylee/0.GWAS"
impute_dir="${Dir}/1.Input/${Sample}/0.Impute"
QC_dir="${Dir}/1.Input/${Sample}/1.QC"
inversion="${Dir}/REFERENCE/inversion.txt"
Rcode="/ichrogene/project/temp/gylee/0.GWAS/Code/R"
p2_Dir="/ichrogene/project/temp/gylee/Singularity/plink2"
Output_dir="${Dir}/2.plink_result/${Sample}"
#=================================================================
ana_dir="${Dir}/2.plink_result/NoImputed_${Sample}"
SAGE_CD="/ichrogene/project/temp/gylee/0.GWAS/Code/SAIGE_code_folder"


echo "=========================="
echo "Start Step0.createSparseGRM.sh"
echo "=========================="
mkdir -p ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LOG
mkdir -p ${Dir}/2.SAIGE_result/${Sample}/result
Log="${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LOG"


#Made by prune data
Rscript ${SAGE_CD}/createSparseGRM.R \
     --plinkFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC/LD_pruned_QC_${Sample} \
     --nThreads=32  \
     --outputPrefix=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/sparseGRM \
     --numRandomMarkerforSparseKin=5000 \
     --relatednessCutoff=0.125 > ${Log}/Stop0.log
