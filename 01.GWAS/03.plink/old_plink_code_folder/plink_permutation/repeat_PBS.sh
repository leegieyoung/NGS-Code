#!/bin/sh

for A in $(seq 1 10000)
do
qsub PBS/PBS_$A.sh
done
