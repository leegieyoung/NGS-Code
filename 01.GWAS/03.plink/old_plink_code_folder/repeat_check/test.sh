#!/bin/sh
Start=$1
number=43
End=$((Start+number))
echo $End
for A in $(seq 1 60)
do
cp PBS_.sh PBS_$A.sh
done



