#!/bin/sh
for A in $(seq 1 20)
do
qsub PBS_repeat_output$A.sh
done
