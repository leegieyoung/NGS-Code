#!/bin/sh
#PBS -V
#PBS -q normal
#PBS -A etc
#PBS -l select=1:ncpus=68:mpiprocs=1:ompthreads=68
#PBS -l walltime=48:00:00
#PBS -o /scratch/hpc46a05/PBS/PBS.OU
#PBS -e /scratch/hpc46a05/PBS/Error
#PBS -N divide13

Code_path="/scratch/hpc46a05/GWAS/Code/plink_code_folder"
module purge
module load python/3.7
. /apps/applications/PYTHON/3.7/etc/profile.d/conda.sh
conda activate plink

#sh ${Code_path}/1-1_AfterBeagle_QC.sh divide13
sh ${Code_path}/1-2_Merge_dataset_kao_2Group.sh snpQC_CD-CODA_divide13 snpQC-CD CODA_QCsex_divide13
sh ${Code_path}/2_Single_CaseControl_analysis_major_kao_QC.sh Imputed_snpQC_CD-CODA_divide13
