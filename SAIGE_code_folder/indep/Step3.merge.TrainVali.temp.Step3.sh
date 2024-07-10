#!/bin/sh
if [ $# -ne 3 ];then
        echo "Please enter Sample_Name, thread(s)"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
#beagle 진행시Same position 제거를 안했을 경우 사용하는 코드
Sample=$1
Thread=$2
TrainVali=$3
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
SAIGE_dir="${Dir}/2.SAIGE_result/${TrainVali}/${Sample}"

mkdir -p ${SAIGE_dir}/sparseGRM/LOG
mkdir -p ${SAIGE_dir}/result
Log="${SAIGE_dir}/sparseGRM/LOG/"
export RESULT_DIR="${SAIGE_dir}/result/"
#echo ${ANA_DIR}
export SAMPLE="${Sample}"
Rscript ${Dir}/Code/R/Manhattan_plot.SAIGE.R



