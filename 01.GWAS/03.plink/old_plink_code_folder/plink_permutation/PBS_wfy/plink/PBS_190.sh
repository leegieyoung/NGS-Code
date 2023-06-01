#!/bin/sh
#PBS -q normal
#PBS -A etc
#PBS -l select=1:ncpus=68:mpiprocs=1:ompthreads=68
#PBS -l walltime=48:00:00
#PBS -o /scratch/hpc46a05/PBS.OU
#PBS -e /scratch/hpc46a05/Error
#PBS -N 190

module purge
module load python/3.7
. /apps/applications/PYTHON/3.7/etc/profile.d/conda.sh
conda activate plink


Sample_path="/scratch/hpc46a05/GWAS/result/AD/maf005/Imputed_BSAD_199_analysis_folder/merge"
Result_folder="/scratch/hpc46a05/GWAS/result/AD/maf005/Imputed_BSAD_199_analysis_folder/merge/permutation/result"
SummRes_folder="/scratch/hpc46a05/GWAS/result/AD/maf005/Imputed_BSAD_199_analysis_folder/merge/permutation/summary_result"
INPUT=190
one=1
multi=$((INPUT-one))
num=50
before_Start=$((multi*num))
Start=$((before_Start+one))
End=$((INPUT*num))

for A in $(seq $Start $End)
do
plink \
 --bed ${Sample_path}/Imputed_BSAD_199_NoNA.bed \
 --bim ${Sample_path}/Imputed_BSAD_199_NoNA.bim \
 --fam /scratch/hpc46a05/GWAS/result/AD/maf005/Imputed_BSAD_199_analysis_folder/merge/permutation/Sample/$A.fam \
 --ci 0.95 \
 --covar ${Sample_path}/logistic.covar_mds_age.txt \
 --keep-allele-order \
 --logistic \
 --hide-covar \
 --sex \
 --out ${Result_folder}/check${A}

awk '$12 < 1e-3 {print $2, $12}' ${Result_folder}/check${A}.assoc.logistic > ${SummRes_folder}/1e3_${A}.txt
sed -i '1i\SNP P' ${SummRes_folder}/1e3_${A}.txt
done