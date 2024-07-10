#!/bin/sh
if [ $# -ne 4 ];then
        echo "Please enter Sample_Name"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
#beagle 진행시Same position 제거를 안했을 경우 사용하는 코드
Sample=$1
Thread=$2
MEM=$3
Cohort=$4
Dir="/ichrogene/project/temp/gylee/0.GWAS"
p2_Dir="/ichrogene/project/temp/gylee/Singularity/plink2"
Inversion="/ichrogene/project/temp/gylee/0.GWAS/REFERENCE/inversion.txt"
Output_dir="${Dir}/2.plink_result/${Sample}_${Cohort}"
Singularity_Dir="/ichrogene/project/temp/gylee/Singularity"
data_Dir="/ichrogene/project/temp/gylee/1.WES/REFERENCE/snpEff/data"
REFERENCE="/ichrogene/project/temp/gylee/1.WES/REFERENCE/"

mkdir -p ${Output_dir} 

export RESULT_DIR="${Output_dir}/"
#echo ${ANA_DIR}
export SAMPLE="${Sample}"

awk '{if($16 < 0.01) print $3}' ${Output_dir}/gachon.Pheno.glm.logistic.hybrid > ${Output_dir}/1e-2.snp.list

${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe_rmoutlier \
 --extract ${Output_dir}/1e-2.snp.list \
 --recode vcf bgz \
 --out ${Output_dir}/1e-2.snp

singularity exec --bind /ichrogene/:/ichrogene/ ${Singularity_Dir}/PRS.sif bcftools \
 index -f --threads ${Thread} \
 ${Output_dir}/1e-2.snp.vcf.gz
#snpEff
java -Xms2g -Xmx2g -jar ${Singularity_Dir}/snpEff/snpEff.jar eff \
 -dataDir ${data_Dir} -v GRCh38.86 \
 ${Output_dir}/1e-2.snp.vcf.gz > ${Output_dir}/1e-2.snp.snpEff.vcf
singularity exec --bind /ichrogene/:/ichrogene/ ${Singularity_Dir}/PRS.sif bcftools \
 view --threads ${Thread} ${Output_dir}/1e-2.snp.snpEff.vcf \
 -Oz -o ${Output_dir}/1e-2.snp.snpEff.vcf.gz
singularity exec --bind /ichrogene/:/ichrogene/ ${Singularity_Dir}/PRS.sif bcftools \
 index -f --threads ${Thread} ${Output_dir}/1e-2.snp.snpEff.vcf.gz

#SnpSift
java -Xms2g -Xmx2g -jar ${Singularity_Dir}/snpEff/SnpSift.jar annotate \
 ${REFERENCE}/dbsnp/dbsnp151.GRCh38.p7.chrM_edit.vcf.gz \
 ${Output_dir}/1e-2.snp.snpEff.vcf.gz > ${Output_dir}/1e-2.snp.snpEff.dbsnp.vcf

grep -v '^##' ${Output_dir}/1e-2.snp.snpEff.dbsnp.vcf > ${Output_dir}/temp.vcf
singularity exec --bind /ichrogene/:/ichrogene/ ${Singularity_Dir}/GWAS.sif Rscript /ichrogene/project/temp/gylee/Code/snpeff_code_folder/Find.Gene.R
rm ${Output_dir}/temp.vcf

awk '{if($16 < 0.01) print $0}' ${Output_dir}/gachon.Pheno.glm.logistic.hybrid > ${Output_dir}/temp.glm
head -n 1 ${Output_dir}/gachon.Pheno.glm.logistic.hybrid > ${Output_dir}/head
cat ${Output_dir}/head ${Output_dir}/temp.glm > ${Output_dir}/1e-2.snp.glm
rm ${Output_dir}/temp.glm


