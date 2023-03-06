#!/bin/sh
for A in $(cat qsub.list)
do
qsub $A
done
