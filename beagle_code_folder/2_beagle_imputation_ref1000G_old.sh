#!/bin/bash
#raw 데이터는 모두 준비가 되어 있어야 함
#raw.vcf 면 됨
#Sample_folder수정필요#############필독#########
if [ $# -ne 2 ];then
	echo "Please enter ##"
		exit
fi
input=$1
Sample=${input}
Chrom=$2
Beagle_Code_folder="/mnt/data/gylee/0.GWAS/Code/beagle_code_folder"
#Sample_folder="/scratch/hpc46a05/GWAS/Eagle2_result/2.${Sample}"
Sample_folder="/mnt/data/gylee/0.GWAS/2.beagle_result/1.${Sample}"
Reference_folder="/mnt/data/gylee/0.GWAS/beagle5.2/1000_Genomes_phase3_v5a/b37.vcf"
Reference_bref3="/mnt/data/gylee/0.GWAS/beagle5.2/1000_Genomes_phase3_v5a/b37.bref3"
Result_folder="/mnt/data/gylee/0.GWAS/2.beagle_result"


mkdir ${Result_folder}/2.${input}_conform-gt_${Sample}
/usr/lib/jdk-16.0.2/bin/java -Xmx64g -jar ${Beagle_Code_folder}/conform-gt.24May16.cee.jar \
 ref=${Reference_folder}/chr${Chrom}.1kg.phase3.v5a.vcf.gz \
 gt=${Sample_folder}/CHR${Chrom}_${Sample}_noIndel_rmDupID.vcf.gz \
 chrom=${Chrom} \
 out=${Result_folder}/2.${input}_conform-gt_${Sample}/raw_CHR${Chrom}_${Sample}_noIndel_rmDupID_beagle \
 match=POS

mkdir ${Result_folder}/3.${input}_imputed_${Sample}
/usr/lib/jdk-16.0.2/bin/java -Xmx64g -jar ${Beagle_Code_folder}/beagle.18May20.d20.jar \
 ref=${Reference_bref3}/chr${Chrom}.1kg.phase3.v5a.b37.bref3 \
 gt=${Result_folder}/2.${input}_conform-gt_${Sample}/raw_CHR${Chrom}_${Sample}_noIndel_rmDupID_beagle.vcf.gz \
 chrom=${Chrom} \
 out=${Result_folder}/3.${input}_imputed_${Sample}/CHR${Chrom}_${Sample}_noIndel_rmDupID_beagle \
 impute=true ap=true window=30

#/usr/lib/jdk-16.0.2/bin/java -Xmx64g -jar ${Beagle_Code_folder}/beagle.22Jul22.46e.jar \
