#!/bin/sh

for A in $(seq 1 8)
do
one=1
num=232
Start=$((${num}*(${A}-${one})+${one}))
End=$((${num}*${A}))

echo ${Start}
echo ${End}

#cp EWparsing(1,250)

cp PGS_EW_step_base.R PGS_EW_step${A}.R
sed -i "s/EWparsing(Start,End)/EWparsing(${Start},${End})/g" PGS_EW_step${A}.R 
done
