#!/bin/bash
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
p2_Dir="/mnt/nas/gylee/Singurality/plink2"
Pcode="${Dir}/Code/plink_code_folder"
Scode="${Dir}/Code/SAIGE_code_folder"
#=================================================================
 
mkdir -p ${Dir}/2.plink_result/NoImputed_${Sample}/temp
mkdir -p ${Dir}/2.plink_result/NoImputed_${Sample}/result
ana_dir="${Dir}/2.plink_result/NoImputed_${Sample}"

singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif bash ${Pcode}/2_preSAIGE.sh ${Sample} ${Thread}
bash ${Pcode}/3_plink2_CaseControl_analysis_major_kao_QC_OR_NoImpute_Memory.sh ${Sample} ${Thread}
singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/SAIGE/SAIGE_docker_R36_1.1.9.sif bash ${Scode}/Step0.createSparseGRM.sh ${Sample} ${Thread}
#singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/SAIGE/SAIGE_docker_R36_1.1.9.sif bash ${Scode}/Step1.fitNULLGLMM.sh ${Sample} ${Thread}
#singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/SAIGE/SAIGE_docker_R36_1.1.9.sif bash ${Scode}/Step2.SPAtests.sh ${Sample} ${Thread}
