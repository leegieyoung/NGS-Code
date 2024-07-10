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
 
mkdir ${SAIGE_dir}/result/temp/
for A in $(seq 2 22)
do
sed '1,1d' ${SAIGE_dir}/result/chr${A}.${Sample}.LDpruned.total.SAIGE.test.markers.txt > ${SAIGE_dir}/result/temp/chr${A}.merge
done

cat ${SAIGE_dir}/result/chr1.${Sample}.LDpruned.total.SAIGE.test.markers.txt \
 ${SAIGE_dir}/result/temp/chr2.merge \
 ${SAIGE_dir}/result/temp/chr3.merge \
 ${SAIGE_dir}/result/temp/chr4.merge \
 ${SAIGE_dir}/result/temp/chr5.merge \
 ${SAIGE_dir}/result/temp/chr6.merge \
 ${SAIGE_dir}/result/temp/chr7.merge \
 ${SAIGE_dir}/result/temp/chr8.merge \
 ${SAIGE_dir}/result/temp/chr9.merge \
 ${SAIGE_dir}/result/temp/chr10.merge \
 ${SAIGE_dir}/result/temp/chr11.merge \
 ${SAIGE_dir}/result/temp/chr12.merge \
 ${SAIGE_dir}/result/temp/chr13.merge \
 ${SAIGE_dir}/result/temp/chr14.merge \
 ${SAIGE_dir}/result/temp/chr15.merge \
 ${SAIGE_dir}/result/temp/chr16.merge \
 ${SAIGE_dir}/result/temp/chr17.merge \
 ${SAIGE_dir}/result/temp/chr18.merge \
 ${SAIGE_dir}/result/temp/chr19.merge \
 ${SAIGE_dir}/result/temp/chr20.merge \
 ${SAIGE_dir}/result/temp/chr21.merge \
 ${SAIGE_dir}/result/temp/chr22.merge > ${SAIGE_dir}/result/temp/merge.result
awk '{print $1,$2,$3,$4,$5,$7,$8,$9,$10,$13,$15,$16,$17,$18,$19}' ${SAIGE_dir}/result/temp/merge.result > ${SAIGE_dir}/result/merge.result

rm -rf ${SAIGE_dir}/result/temp/

export RESULT_DIR="${SAIGE_dir}/result/"
#echo ${ANA_DIR}
export SAMPLE="${Sample}"
Rscript ${Dir}/Code/R/Manhattan_plot.SAIGE.R



