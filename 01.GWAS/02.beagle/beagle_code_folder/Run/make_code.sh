#!/bin/sh
for A in $(cat qsub.list)
do
sed -i 's/divide2/divide1/g' $A 
#sed -i 's/\/scratch\/hpc46a05\/PBS.OU/\/scratch\/hpc46a05\/PBS\/PBS.OU/g' $A
#sed -i 's/\/scratch\/hpc46a05\/Error/\/scratch\/hpc46a05\/PBS\/Error/g' $A
done
