#!/bin/sh
#for A in $(seq 1 60)
#do
#qsub PBS_$A.sh
#done

for B in $(seq 1 7)
do
qsub PBS_flat$B.sh
done
