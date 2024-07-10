#!/bin/sh
for A in $(seq 1 22)
do
sh weCall.sh ${A} &
done
