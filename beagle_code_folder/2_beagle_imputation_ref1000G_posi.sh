#!/bin/bash
#raw 데이터는 모두 준비가 되어 있어야 함
#raw.vcf 면 됨
#Sample_folder수정필요#############필독#########
if [ $# -ne 4 ];then
	echo "Please enter input_name, chr, start and end."
		exit
fi
input=$1
chr=$2
Start=$3
End=$4
Beagle_Code_folder="/mnt/nas/gylee/0.GWAS/Code/beagle_code_folder"
Sample_folder="/mnt/nas/gylee/0.GWAS/2.beagle_result/1.${input}"
Reference_folder="/mnt/nas/gylee/0.GWAS/REFERENCE/1000_Genomes_phase3_v5a/b37.vcf"
Reference_bref3="/mnt/nas/gylee/0.GWAS/REFERENCE/1000_Genomes_phase3_v5a/b37.bref3"
Result_folder="/mnt/nas/gylee/0.GWAS/2.beagle_result"


mkdir ${Result_folder}/2.conform-gt_${input}
/usr/lib/jdk-16.0.2/bin/java -Xmx64g -jar ${Beagle_Code_folder}/conform-gt.24May16.cee.jar \
 ref=${Reference_folder}/chr${chr}.1kg.phase3.v5a.vcf.gz \
 gt=${Sample_folder}/CHR${chr}_${input}_noIndel_rmDupID.vcf.gz \
 chrom=${chr} \
 out=${Result_folder}/2.conform-gt_${input}/raw_CHR${chr}_${input}_noIndel_rmDupID_beagle \
 match=POS

mkdir ${Result_folder}/3.imputed_${input}
/usr/lib/jdk-16.0.2/bin/java -Xmx64g -jar ${Beagle_Code_folder}/beagle.18May20.d20.jar \
 ref=${Reference_bref3}/chr${chr}.1kg.phase3.v5a.b37.bref3 \
 gt=${Result_folder}/2.${input}_conform-gt_${input}/raw_CHR${chr}_${input}_noIndel_rmDupID_beagle.vcf.gz \
 chrom=${chr}:${Start}-${End} \
 out=${Result_folder}/3.imputed_${input}/CHR${chr}_${input}_noIndel_rmDupID_beagle \
 impute=true ap=true window=30

#/usr/lib/jdk-16.0.2/bin/java -Xmx64g -jar ${Beagle_Code_folder}/beagle.22Jul22.46e.jar \
