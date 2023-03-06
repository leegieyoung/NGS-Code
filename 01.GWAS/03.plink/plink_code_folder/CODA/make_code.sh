#!/bin/sh
for A in $(seq 1 18)
do
cp PBS_merge_CD-CODA_QCsex_divide1.sh PBS_merge_CD-CODA_QCsex_divide${A}.sh
sed -i "s/divide1/divide${A}/g" PBS_merge_CD-CODA_QCsex_divide${A}.sh
done
