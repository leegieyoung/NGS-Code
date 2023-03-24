#!/bin/sh
#PBS -V
#PBS -q normal
#PBS -A etc
#PBS -l select=1:ncpus=68:mpiprocs=1:ompthreads=68
#PBS -l walltime=48:00:00
#PBS -o /scratch/x1997a11/PBS.OU
#PBS -e /scratch/x1997a11/Error
#PBS -N 6

Code_path="/scratch/x1997a11/GWAS/pdxen_AD/Code_folder/keep_allele_order_version"
module purge
module load python/3.7
. /apps/applications/PYTHON/3.7/etc/profile.d/conda.sh
conda activate plink

sh ${Code_path}/1-1_Imputed_vcf_to_bfile.sh AD_14 6
sh ${Code_path}/1-1_Imputed_vcf_to_bfile.sh AD_185 6
