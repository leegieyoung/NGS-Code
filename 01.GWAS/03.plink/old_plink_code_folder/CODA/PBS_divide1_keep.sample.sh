#!/bin/sh
#PBS -V
#PBS -q normal
#PBS -A etc
#PBS -l select=1:ncpus=68:mpiprocs=1:ompthreads=68
#PBS -l walltime=48:00:00
#PBS -o /scratch/hpc46a05/PBS/PBS.OU
#PBS -e /scratch/hpc46a05/PBS/Error
#PBS -N div_CODA

#!/bin/sh
GWAS_path="/scratch/hpc46a05/GWAS/result"
code_path="/scratch/hpc46a05/GWAS/Code/plink_code_folder"
inversion="/scratch/hpc46a05/REFERENCE/inversion.txt"
Rcode="/scratch/hpc46a05/GWAS/Rcode"
#=================================================================

module purge
module load python/3.7
. /apps/applications/PYTHON/3.7/etc/profile.d/conda.sh
conda activate plink


for A in  $(seq 1 18)
do
mkdir -p /scratch/hpc46a05/GWAS/result/QC_Imputed_CODA_QCsex_divide${A}/MaMi
plink --bfile /scratch/hpc46a05/GWAS/result/QC_Imputed_CODA_QCsex/MaMi/QC_Imputed_CODA_QCsex \
 --keep-allele-order \
 --indiv-sort n \
 --make-bed \
 --keep /scratch/hpc46a05/GWAS/result/CODA/divide${A}.sample \
 --out /scratch/hpc46a05/GWAS/result/QC_Imputed_CODA_QCsex_divide${A}/MaMi/QC_Imputed_CODA_QCsex_divide${A}
done
