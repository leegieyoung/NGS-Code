##!/bin/sh
if [ $# -ne 3 ];then
        echo "Please enter Sample_Name"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
#beagle 진행시Same position 제거를 안했을 경우 사용하는 코드
Sample=$1
Thread=$2
MEM=$3
Dir="/mnt/nas/gylee/0.GWAS"
QC_dir="${Dir}/1.Input/${Sample}/1.QC"
MaMi_dir="${Dir}/1.Input/${Sample}/2.MaMi"
p2_Dir="/mnt/nas/gylee/Singurality/plink2"
Inversion="/mnt/nas/gylee/0.GWAS/REFERENCE/inversion.txt"
Output_dir="${Dir}/2.plink_result/${Sample}"

mkdir -p ${QC_dir}
mkdir -p ${Output_dir} 

export QC_DIR="${QC_dir}/"
export RESULT_DIR="${Output_dir}/"
echo ${ANA_DIR}
export SAMPLE="${Sample}"


singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/SAIGE/SAIGE_docker_R36.sif bash /mnt/nas/gylee/0.GWAS/Code/SAIGE_code_folder/Step0.createSparseGRM.sh ${Sample} ${Thread}
singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/SAIGE/SAIGE_docker_R36.sif bash /mnt/nas/gylee/0.GWAS/Code/SAIGE_code_folder/Step1.fitNULLGLMM.sh ${Sample} ${Thread}
singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/SAIGE/SAIGE_docker_R36.sif bash /mnt/nas/gylee/0.GWAS/Code/SAIGE_code_folder/Step2.SPAtests_exceptchr2.sh ${Sample} ${Thread}
singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/SAIGE/SAIGE_docker_R36.sif bash /mnt/nas/gylee/0.GWAS/Code/SAIGE_code_folder/Step2.SPAtests_chr2.sh ${Sample} ${Thread}

