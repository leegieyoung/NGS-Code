#!/bin/sh
for A in $(seq 1 15)
do
qsub PBS_flat$A.sh
done
