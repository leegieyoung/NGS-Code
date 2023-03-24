#!/bin/sh
for A in $(cat /scratch/x1997a11/GWAS/pdxen_AD/Code_folder/keep_allele_order_version/PBS_BSAD/repeat.list)
do
qsub ${A}
done
