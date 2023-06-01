#!/bin/sh

for A in $(cat qdel.list)
do
qdel $A
done
