#!/bin/sh
#PBS -V
#PBS -q normal
#PBS -A etc
#PBS -l select=1:ncpus=68:mpiprocs=1:ompthreads=68
#PBS -l walltime=48:00:00
#PBS -o /scratch/x1997a11/PBS.OU
#PBS -e /scratch/x1997a11/Error
#PBS -N 9001-9500
for A in $(seq 9001 9500)
do
Path="/scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check/plink_pphe"
echo "check_${A}" > ${Path}/check$A/${A}.txt

grep '1:247601357'  ${Path}/check$A/check${A}.assoc.logistic > ${Path}/check$A/check${A}.rs12143966
paste ${Path}/check$A/${A}.txt ${Path}/check$A/check${A}.rs12143966 >> /scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check/rs12143966
grep '1:247603463'  ${Path}/check$A/check${A}.assoc.logistic > ${Path}/check$A/check${A}.rs4925659
paste ${Path}/check$A/${A}.txt ${Path}/check$A/check${A}.rs4925659 >> /scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check/rs4925659
grep '1:247604447'  ${Path}/check$A/check${A}.assoc.logistic > ${Path}/check$A/check${A}.rs111307268
paste ${Path}/check$A/${A}.txt ${Path}/check$A/check${A}.rs111307268 >> /scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check/rs111307268
grep '1:247607052'  ${Path}/check$A/check${A}.assoc.logistic > ${Path}/check$A/check${A}.rs10159239
paste ${Path}/check$A/${A}.txt ${Path}/check$A/check${A}.rs10159239 >> /scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check/rs10159239
grep '1:247607642'  ${Path}/check$A/check${A}.assoc.logistic > ${Path}/check$A/check${A}.rs12130711
paste ${Path}/check$A/${A}.txt ${Path}/check$A/check${A}.rs12130711 >> /scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check/rs12130711
grep '1:247610411'  ${Path}/check$A/check${A}.assoc.logistic > ${Path}/check$A/check${A}.rs56383829
paste ${Path}/check$A/${A}.txt ${Path}/check$A/check${A}.rs56383829 >> /scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check/rs56383829
grep '1:247612036'  ${Path}/check$A/check${A}.assoc.logistic > ${Path}/check$A/check${A}.rs10754558
paste ${Path}/check$A/${A}.txt ${Path}/check$A/check${A}.rs10754558 >> /scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check/rs10754558
grep '1:247612295'  ${Path}/check$A/check${A}.assoc.logistic > ${Path}/check$A/check${A}.rs10802502
paste ${Path}/check$A/${A}.txt ${Path}/check$A/check${A}.rs10802502 >> /scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check/rs10802502
grep '1:247612435'  ${Path}/check$A/check${A}.assoc.logistic > ${Path}/check$A/check${A}.rs4925547
paste ${Path}/check$A/${A}.txt ${Path}/check$A/check${A}.rs4925547 >> /scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check/rs4925547
grep '1:247612562'  ${Path}/check$A/check${A}.assoc.logistic > ${Path}/check$A/check${A}.rs10925027
paste ${Path}/check$A/${A}.txt ${Path}/check$A/check${A}.rs10925027 >> /scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check/rs10925027
grep '1:247613321'  ${Path}/check$A/check${A}.assoc.logistic > ${Path}/check$A/check${A}.rs10754559
paste ${Path}/check$A/${A}.txt ${Path}/check$A/check${A}.rs10754559 >> /scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check/rs10754559
grep '1:247613348'  ${Path}/check$A/check${A}.assoc.logistic > ${Path}/check$A/check${A}.rs10733111
paste ${Path}/check$A/${A}.txt ${Path}/check$A/check${A}.rs10733111 >> /scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check/rs10733111
grep '1:247613436'  ${Path}/check$A/check${A}.assoc.logistic > ${Path}/check$A/check${A}.rs10802503
paste ${Path}/check$A/${A}.txt ${Path}/check$A/check${A}.rs10802503 >> /scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check/rs10802503
grep '1:247613888'  ${Path}/check$A/check${A}.assoc.logistic > ${Path}/check$A/check${A}.rs10733112
paste ${Path}/check$A/${A}.txt ${Path}/check$A/check${A}.rs10733112 >> /scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check/rs10733112
grep '1:247614617'  ${Path}/check$A/check${A}.assoc.logistic > ${Path}/check$A/check${A}.rs4925663
paste ${Path}/check$A/${A}.txt ${Path}/check$A/check${A}.rs4925663 >> /scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check/rs4925663
grep '1:247614896'  ${Path}/check$A/check${A}.assoc.logistic > ${Path}/check$A/check${A}.rs11583410
paste ${Path}/check$A/${A}.txt ${Path}/check$A/check${A}.rs11583410 >> /scratch/x1997a11/GWAS/pdxen_AD/result_folder/repeat_check/rs11583410
done

