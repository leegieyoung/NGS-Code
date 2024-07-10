#!/bin/sh
#PBS -V
#PBS -q normal
#PBS -A etc
#PBS -l select=1:ncpus=68:mpiprocs=1:ompthreads=68
#PBS -l walltime=48:00:00
#PBS -o /scratch/hpc46a05/PBS/PBS.OU
#PBS -e /scratch/hpc46a05/PBS/Error
#PBS -N CD-CODA-merge

Code_path="/scratch/hpc46a05/GWAS/Code/plink_code_folder"
module purge
module load python/3.7
. /apps/applications/PYTHON/3.7/etc/profile.d/conda.sh
conda activate plink

#sh ${Code_path}/1-2_Merge_dataset_kao_2Group.sh snpQC_pcaQC-CD-CODA snpQC-CD snpQC-800CODA
sh ${Code_path}/1-2_Merge_dataset_kao_2Group.sh snpQC_pcaQC-CD-33702CODA snpQC-CD CODA_QCsex
