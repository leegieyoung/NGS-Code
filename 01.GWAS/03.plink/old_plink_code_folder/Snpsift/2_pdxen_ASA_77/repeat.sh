#!/bin/sh
for A in $(cat repeat.list)
do
qsub ${A}
done
