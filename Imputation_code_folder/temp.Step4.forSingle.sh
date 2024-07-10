#!/bin/sh
if [ $# -ne 2 ];then
        echo "Please enter INPUT file name & Chromosme(Only number)"
               exit
fi

Code_Dir="/ichrogene/project/temp/gylee/0.GWAS/Code/Imputation_code_folder"
Dir="/ichrogene/project/temp/gylee/0.GWAS/05.Imputation" #Change to your path
Image="/ichrogene/project/temp/gylee/Singularity"  #Change to your image file path
p2_Dir="/ichrogene/project/temp/gylee/Singularity/plink2"
scratch="/ichrogene"
ref_Dir="/ichrogene/project/temp/gylee/0.GWAS/REFERENCE/Imputation/reference"

Input=$1
Chr=$2
Thread=3
Chr5=5
Chr6=6
mkdir ${Code_Dir}/LOG/

minimac_start() {
 local Input=$1
 local Chr=$2
 local Thread=$3
	bcftools index -f ${Dir}/1.shapeit/imputation_${Input}_${Chr}_shapeit.vcf.gz ;

	minimac4 --format DS,GT,GP \
        	--referenceEstimates \
	        --mapFile ${ref_Dir}/map/chr${Chr}.b37.gmap.gz \
	        --refHaps ${ref_Dir}/minimac/${Chr}.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz \
	        --haps ${Dir}/1.shapeit/imputation_${Input}_${Chr}_shapeit.vcf.gz \
	        --prefix ${Dir}/2.minimac4/minimac4.1.0/imputation_${Input}_${Chr}_shapeit_minimac4.1.0 --log --cpus ${Thread} 
}

minimac_start ${Input} ${Chr} 1
