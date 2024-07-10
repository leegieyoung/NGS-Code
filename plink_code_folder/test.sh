#!/bin/sh
if [ $# -ne 6 ];then
        echo "Please enter Sample_Name"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
#beagle 진행시Same position 제거를 안했을 경우 사용하는 코드
Sample=$1
Thread=$2
MEM=$3
CUTOFF=$4
TRAIT=$5
TrainVali=$6
Dir="/ichrogene/project/temp/gylee/0.GWAS"
QC_dir="${Dir}/1.Input/${Sample}/1.QC"
MaMi_dir="${Dir}/1.Input/${Sample}/2.MaMi"
p2_Dir="/ichrogene/project/temp/gylee/Singularity/plink2"
Inversion="/ichrogene/project/temp/gylee/0.GWAS/REFERENCE/inversion.txt"
Output_dir="${Dir}/2.plink_result/${TrainVali}/${Sample}"
SAIGE_dir="${Dir}/2.SAIGE_result/${TrainVali}/${Sample}/result"
mkdir -p ${QC_dir}
mkdir -p ${Output_dir} 

export QC_DIR="${QC_dir}/"
export RESULT_DIR="${Output_dir}/"
echo ${ANA_DIR}
export SAMPLE="${Sample}"
export SAIGE_DIR="${SAIGE_dir}/"
export TRAIT="${TRAIT}"

export CUTOFF=${CUTOFF}

#SHAP
#make CT.vcf
#${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe \
# --extract ${Output_dir}/prune/iMAC.PRSice.p1.r10.window250kb.prune.in \
# --recode vcf \
# --out ${Output_dir}/prune/${TRAIT}.CT

#singularity exec --bind /ichrogene/:/ichrogene/ /ichrogene/project/temp/gylee/Singularity/AI.sif Rscript /ichrogene/project/temp/gylee/0.GWAS/4.PRS/shap.R 


CUTOFF="${TRAIT}.CT.SHAP"
export CUTOFF=${CUTOFF}

Rscript /ichrogene/project/temp/gylee/0.GWAS/Code/R/make_marker_forSAIGE.R
