!/bin/sh
if [ $# -ne 3 ];then
        echo "Please enter Sample_Name"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
#beagle 진행시Same position 제거를 안했을 경우 사용하는 코드
Sample=$1
Thread=$2
MEM=$3
Dir="/mnt/nas/gylee/0.GWAS"
QC_dir="${Dir}/1.Input/${Sample}/1.QC"
MaMi_dir="${Dir}/1.Input/${Sample}/2.MaMi"
p2_Dir="/mnt/nas/gylee/Singurality/plink2"
Inversion="/mnt/nas/gylee/0.GWAS/REFERENCE/inversion.txt"
Output_dir="${Dir}/2.plink_result/${Sample}"

mkdir -p ${QC_dir}
mkdir -p ${Output_dir} 

export QC_DIR="${QC_dir}/"
export RESULT_DIR="${Output_dir}/"
echo ${ANA_DIR}
export SAMPLE="${Sample}"

awk '{print $3}' ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid.FDR > ${Output_dir}/${Sample}.FDR.snplist
awk '{print $3}' ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid.suggest > ${Output_dir}/${Sample}.suggest.snplist
awk '{print $3}' ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid.1e-4 > ${Output_dir}/${Sample}.1e-4.snplist

#FDR.min0.2
mkdir ${Output_dir}/prune
${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe \
 --extract ${Output_dir}/${Sample}.FDR.snplist \
 --snps-only just-acgt \
 --make-pgen \
 --out ${Output_dir}/prune/FDR.just-acgt

${p2_Dir}/plink2 --pfile ${Output_dir}/prune/FDR.just-acgt \
 --indep-pairwise 1000kb 0.2 \
 --pheno ${Output_dir}/${Sample}_pheno.txt \
 --memory ${MEM} \
 --out ${Output_dir}/prune/FDR.min0.2

for A in $(cat ${Output_dir}/prune/FDR.min0.2.prune.in)
do
grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid >> ${Output_dir}/prune/raw.FDR.min0.2.glm
done

head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
cat ${Output_dir}/header ${Output_dir}/prune/raw.FDR.min0.2.glm > ${Output_dir}/prune/FDR.min0.2.glm
rm -rf ${Output_dir}/prune/raw.FDR.min0.2.glm

export CUTOFF="FDR.min0.2"
CUTOFF="FDR.min0.2"
singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker.R
bash /mnt/nas/gylee/0.GWAS/4.PRS/KoGES/calculate_PRS.230713.KOR.sh ${Sample} type2diabet obesity ${CUTOFF}

#FDR.min0.02
${p2_Dir}/plink2 --pfile ${Output_dir}/prune/FDR.just-acgt \
 --indep-pairwise 1000kb 0.02 \
 --pheno ${Output_dir}/${Sample}_pheno.txt \
 --memory ${MEM} \
 --out ${Output_dir}/prune/FDR.min0.02

for A in $(cat ${Output_dir}/prune/FDR.min0.02.prune.in)
do
grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid >> ${Output_dir}/prune/raw.FDR.min0.02.glm
done

head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
cat ${Output_dir}/header ${Output_dir}/prune/raw.FDR.min0.02.glm > ${Output_dir}/prune/FDR.min0.02.glm
rm -rf ${Output_dir}/prune/raw.FDR.min0.02.glm

export CUTOFF="FDR.min0.02"
CUTOFF="FDR.min0.02"
singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker.R
bash /mnt/nas/gylee/0.GWAS/4.PRS/KoGES/calculate_PRS.230713.KOR.sh ${Sample} type2diabet obesity ${CUTOFF}

rm -rf ${Output_dir}/prune/FDR.just-acgt*

#suggest.min0.2
${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe \
 --extract ${Output_dir}/${Sample}.suggest.snplist \
 --snps-only just-acgt \
 --make-pgen \
 --out ${Output_dir}/prune/suggest.just-acgt

${p2_Dir}/plink2 --pfile ${Output_dir}/prune/suggest.just-acgt \
 --indep-pairwise 1000kb 0.2 \
 --pheno ${Output_dir}/${Sample}_pheno.txt \
 --memory ${MEM} \
 --out ${Output_dir}/prune/suggest.min0.2

for A in $(cat ${Output_dir}/prune/suggest.min0.2.prune.in)
do
grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid >> ${Output_dir}/prune/raw.suggest.min0.2.glm
done

head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
cat ${Output_dir}/header ${Output_dir}/prune/raw.suggest.min0.2.glm > ${Output_dir}/prune/suggest.min0.2.glm
rm -rf ${Output_dir}/prune/raw.suggest.min0.2.glm

export CUTOFF="suggest.min0.2"
CUTOFF="suggest.min0.2"
singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker.R
bash /mnt/nas/gylee/0.GWAS/4.PRS/KoGES/calculate_PRS.230713.KOR.sh ${Sample} type2diabet obesity ${CUTOFF}

#suggest.min.0.02
${p2_Dir}/plink2 --pfile ${Output_dir}/prune/suggest.just-acgt \
 --indep-pairwise 1000kb 0.02 \
 --pheno ${Output_dir}/${Sample}_pheno.txt \
 --memory ${MEM} \
 --out ${Output_dir}/prune/suggest.min0.02

for A in $(cat ${Output_dir}/prune/suggest.min0.02.prune.in)
do
grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid >> ${Output_dir}/prune/raw.suggest.min0.02.glm
done

head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
cat ${Output_dir}/header ${Output_dir}/prune/raw.suggest.min0.02.glm > ${Output_dir}/prune/suggest.min0.02.glm
rm -rf ${Output_dir}/prune/raw.suggest.min0.02.glm

export CUTOFF="suggest.min0.02"
CUTOFF="suggest.min0.02"
singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker.R
bash /mnt/nas/gylee/0.GWAS/4.PRS/KoGES/calculate_PRS.230713.KOR.sh ${Sample} type2diabet obesity ${CUTOFF}

rm -rf ${Output_dir}/prune/suggest.just-acgt*

#1e-4.min.0.02
${p2_Dir}/plink2 --pfile ${Output_dir}/imputation_g_m_maf_hwe \
 --extract ${Output_dir}/${Sample}.1e-4.snplist \
 --snps-only just-acgt \
 --make-pgen \
 --out ${Output_dir}/prune/1e-4.just-acgt

${p2_Dir}/plink2 --pfile ${Output_dir}/prune/1e-4.just-acgt \
 --indep-pairwise 1000kb 0.2 \
 --pheno ${Output_dir}/${Sample}_pheno.txt \
 --memory ${MEM} \
 --out ${Output_dir}/prune/e-4.min0.2

for A in $(cat ${Output_dir}/prune/e-4.min0.2.prune.in)
do
grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid >> ${Output_dir}/prune/raw.e-4.min0.2.glm
done

head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
cat ${Output_dir}/header ${Output_dir}/prune/raw.e-4.min0.2.glm > ${Output_dir}/prune/e-4.min0.2.glm
rm -rf ${Output_dir}/prune/raw.e-4.min0.2.glm

export CUTOFF="e-4.min0.2"
CUTOFF="e-4.min0.2"
singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker.R
bash /mnt/nas/gylee/0.GWAS/4.PRS/KoGES/calculate_PRS.230713.KOR.sh ${Sample} type2diabet obesity ${CUTOFF}


#1e-4.min.0.02
${p2_Dir}/plink2 --pfile ${Output_dir}/prune/1e-4.just-acgt \
 --indep-pairwise 1000kb 0.02 \
 --pheno ${Output_dir}/${Sample}_pheno.txt \
 --memory ${MEM} \
 --out ${Output_dir}/prune/e-4.min0.02

for A in $(cat ${Output_dir}/prune/e-4.min0.02.prune.in)
do
grep -w "${A}" ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid >> ${Output_dir}/prune/raw.e-4.min0.02.glm
done

head -n 1 ${Output_dir}/${Sample}.PHENO1.glm.logistic.hybrid > ${Output_dir}/header
cat ${Output_dir}/header ${Output_dir}/prune/raw.e-4.min0.02.glm > ${Output_dir}/prune/e-4.min0.02.glm
rm -rf ${Output_dir}/prune/raw.e-4.min0.02.glm

export CUTOFF="e-4.min0.02"
CUTOFF="e-4.min0.02"
singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/make_marker.R
bash /mnt/nas/gylee/0.GWAS/4.PRS/KoGES/calculate_PRS.230713.KOR.sh ${Sample} type2diabet obesity ${CUTOFF}
bash /mnt/nas/gylee/0.GWAS/4.PRS/KoGES/calculate_PRS.230713.KOR.sh ${Sample} type2diabet obesity ${CUTOFF} 

rm -rf ${Output_dir}/prune/1e-4.just-acgt*

