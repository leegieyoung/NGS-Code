#!/bin/sh
#PBS -V
#PBS -q long
#PBS -A etc
#PBS -l select=1:ncpus=68:mpiprocs=1:ompthreads=68
#PBS -l walltime=120:00:00
#PBS -o /scratch/x1997a11/PBS.OU
#PBS -e /scratch/x1997a11/Error
#PBS -N 3488
module purge
module load python/3.7
. /apps/applications/PYTHON/3.7/etc/profile.d/conda.sh
conda activate plink

sh /scratch/x1997a11/GWAS/pdxen_AD/Code_folder/repeat_check/longversion_repeat_check_start_end.sh 3488
sh /scratch/x1997a11/GWAS/pdxen_AD/Code_folder/repeat_check/longversion_repeat_check_start_end.sh 3488
