#!/bin/sh
for A in $(seq 1 10)
do
qsub PBS_long$A.sh
done
