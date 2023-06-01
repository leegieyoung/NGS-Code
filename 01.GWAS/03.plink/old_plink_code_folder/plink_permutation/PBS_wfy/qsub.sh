#!/bin/sh
for A in $(cat qsub.list)
do
qsub $A
sleep 1
done
