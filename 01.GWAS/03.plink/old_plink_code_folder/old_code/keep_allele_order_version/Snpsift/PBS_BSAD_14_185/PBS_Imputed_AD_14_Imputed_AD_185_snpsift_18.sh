#!/bin/sh
#PBS -V
#PBS -q normal
#PBS -A etc
#PBS -l select=1:ncpus=68:mpiprocs=1:ompthreads=68
#PBS -l walltime=48:00:00
#PBS -o /scratch/x1997a11/PBS.OU
#PBS -e /scratch/x1997a11/Error
#PBS -N 18
module purge
module load python/3.7
. /apps/applications/PYTHON/3.7/etc/profile.d/conda.sh
conda activate plink
Code_path="/scratch/x1997a11/GWAS/pdxen_AD/Code_folder/keep_allele_order_version"

sh ${Code_path}/2-3_Snpsift.sh Imputed_AD_14_Imputed_AD_185 18