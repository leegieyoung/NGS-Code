#!/bin/sh
#PBS -q normal
#PBS -A etc
#PBS -l select=1:ncpus=68:mpiprocs=1:ompthreads=68
#PBS -l walltime=48:00:00
#PBS -o /scratch/hpc46a05/PBS.OU
#PBS -e /scratch/hpc46a05/Error
#PBS -N 3

module purge
module load python/3.7
. /apps/applications/PYTHON/3.7/etc/profile.d/conda.sh
conda activate plink


Sample_path="/scratch/hpc46a05/GWAS/result/AD/maf005/Imputed_BSAD_199_analysis_folder/merge"
Result_folder="/scratch/hpc46a05/GWAS/result/AD/maf005/Imputed_BSAD_199_analysis_folder/merge/permutation/result"
SummRes_folder="/scratch/hpc46a05/GWAS/result/AD/maf005/Imputed_BSAD_199_analysis_folder/merge/permutation/summary_result"
INPUT=3
one=1
multi=$((INPUT-one))
num=50
before_Start=$((multi*num))
Start=$((before_Start+one))
End=$((INPUT*num))

for A in $(seq $Start $End)
do

awk '{print $2, $12}' ${Result_folder}/check${A}.assoc.logistic > ${SummRes_folder}/${A}.txt
sed -i '1i\SNP P' ${SummRes_folder}/${A}.txt
done