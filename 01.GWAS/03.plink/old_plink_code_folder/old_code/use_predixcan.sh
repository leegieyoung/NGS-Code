#!/bin/bash

tool=/data/yebin/twas/prediXcan/tool
model=/data/yebin/twas/prediXcan/model
weights=("gtex_v7_Colon_Sigmoid_imputed_europeans_tw_0.5_signif.db" "gtex_v7_Colon_Transverse_imputed_europeans_tw_0.5_signif.db" "gtex_v7_Whole_Blood_imputed_europeans_tw_0.5_signif.db")
bfile=/data/yebin/gwas/CD_894.PlEx.imputed.all_ChrPos_unique_noX_appendRS_nodup
pheno=/data/yebin/gwas/UC_twas/CDUC291_X_all_pheno.txt #CD240_pheno.prn
pheno_pre=/data/yebin/gwas/UC_twas/CDUC291_X_all_pheno_pre.txt
pout=/data/yebin/twas/prediXcan/input/CDUC291 #CD_894
bout=/data/yebin/twas/prediXcan/input/CDUC291_cut #CD_240
input=/data/yebin/twas/prediXcan/input/UC #drndr
output=/data/yebin/twas/prediXcan/output/v7
filter=/data/yebin/twas/prediXcan/input/CDUC291_filter.txt
:<<END
#if I have a bed, bim, fam file with more samples than the samples to be used
plink --bfile $bfile --recode --pheno $pheno --allow-no-sex --out $pout #missing pheno = -9
awk '$6 != -9' $pout'.ped' > $bout'.ped'
plink --file $bout --allow-no-sex --make-bed --out $bout

#mkdir $input
END
name=$(echo $input | awk -F "/" '{print $NF}')

#convert plink to dosage file
#python2 $tool/convert_plink_to_dosage.py -b $bout -p plink -o $input/$name
#cp $bout'.fam' $input/$name'_sample.txt'

#name=$(echo $input | awk -F "/" '{print $NF}')
for i in {0..2}; do
	for w in ${weights[i]}; do
		weight=`find -name $w`
		w_name=$(echo $weight | awk -F "/" '{print $NF}'| awk -F "_" '{print $3"_"$4i}')
		echo $w_name

#predicting/imputing expression
		#python2 $tool/PrediXcan.py --predict --dosages $input/ --dosages_prefix $name --samples $name'_sample.txt' --weights $weight --output_prefix $output/$name'_'$w_name

#running association with Phenotype
		python2 $tool/PrediXcan.py --assoc --logistic \
		--pheno $pheno_pre --pheno_name Y \
		--pred_exp $output/$name'_'$w_name'_predicted_expression.txt' \
		--filter $filter 1 \
		--output_prefix $output/$name'_'$w_name

	done
done
