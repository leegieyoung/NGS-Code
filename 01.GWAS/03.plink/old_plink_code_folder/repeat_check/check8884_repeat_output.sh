#!/bin/sh
#PBS -V
#PBS -q normal
#PBS -A etc
#PBS -l select=1:ncpus=68:mpiprocs=1:ompthreads=68
#PBS -l walltime=48:00:00
#PBS -o /scratch/x1997a11/PBS.OU
#PBS -e /scratch/x1997a11/Error
#PBS -N 1-500
A=8884
Path="/scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check/plink_pphe"
echo "check_${A}" > ${Path}/check$A/${A}.txt

grep '1:247614896'  ${Path}/check$A/check${A}.assoc.logistic > ${Path}/check$A/check${A}.rs11583410 
paste ${Path}/check$A/${A}.txt ${Path}/check$A/check${A}.rs11583410 >> /scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check/rs11583410

grep '1:247614617'  ${Path}/check$A/check${A}.assoc.logistic > ${Path}/check$A/check${A}.rs4925663
paste ${Path}/check$A/${A}.txt ${Path}/check$A/check${A}.rs4925663 >> /scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check/rs4925663

grep '20:4705718'  ${Path}/check$A/check${A}.assoc.logistic > ${Path}/check$A/check${A}.rs2245220
paste ${Path}/check$A/${A}.txt ${Path}/check$A/check${A}.rs2245220 >> /scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check/rs2245220

grep '20:4706310'  ${Path}/check$A/check${A}.assoc.logistic > ${Path}/check$A/check${A}.rs2756261
paste ${Path}/check$A/${A}.txt ${Path}/check$A/check${A}.rs2756261 >> /scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check/rs2756261

grep '10:101980355'  ${Path}/check$A/check${A}.assoc.logistic > ${Path}/check$A/check${A}.rs2230803
paste ${Path}/check$A/${A}.txt ${Path}/check$A/check${A}.rs2230803 >> /scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check/rs2230803
