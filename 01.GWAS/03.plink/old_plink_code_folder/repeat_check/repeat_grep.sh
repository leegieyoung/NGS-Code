#!/bin/sh
for A in $(seq 1 20)
do
qsub repeat_output$A.sh
done
