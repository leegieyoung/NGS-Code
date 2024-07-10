#!/bin/sh
if [ $# -ne 2 ];then
        echo "Please enter Sample_Name, thread(s)"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
#beagle 진행시Same position 제거를 안했을 경우 사용하는 코드
Sample=$1
Thread=$2
Dir="/mnt/nas/gylee/0.GWAS"
impute_dir="${Dir}/1.Input/${Sample}/0.Impute"
QC_dir="${Dir}/1.Input/${Sample}/1.QC"
MaMi_dir="${Dir}/1.Input/${Sample}/2.MaMi"
Case_pheno="${ana_dir}/temp/11.Case_pheno.txt"
Control_pheno="${ana_dir}/temp/11.Control_pheno.txt"
inversion="${Dir}/REFERENCE/inversion.txt"
Rcode="/mnt/nas/gylee/0.GWAS/Code/R"
p2_Dir="/mnt/nas/gylee/Singurality/plink2"
#=================================================================
ana_dir="${Dir}/2.plink_result/NoImputed_${Sample}"
SAGE_CD="/mnt/nas/gylee/0.GWAS/Code/SAIGE_code_folder"

echo "=========================="
echo "Start Step1.fitNULLGLMM.sh"
echo "=========================="


mkdir ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LOG
Log="${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LOG/"

for A in $(seq 1 22)
do
Rscript ${SAGE_CD}/step1_fitNULLGLMM.R     \
       --sparseGRMFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/sparseGRM_relatednessCutoff_0.125_5000_randomMarkersUsed.sparseGRM.mtx \
        --sparseGRMSampleIDFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/sparseGRM_relatednessCutoff_0.125_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
        --plinkFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${A}.LD_pruned_QC_${Sample}  \
        --phenoFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/${Sample}.covariate.txt \
       --useSparseGRMtoFitNULL=TRUE \
       --skipVarianceRatioEstimation=FALSE \
        --phenoCol=PHENO \
       --nThreads=${Thread} \
        --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,AGE,SEX \
        --qCovarColList=SEX \
        --sampleIDColinphenoFile=IID \
        --traitType=binary \
        --outputPrefix=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/${Sample}/sparseGRM/fitNGLMM.chr${A}.${Sample}.LDpruned.total \
        --IsOverwriteVarianceRatioFile=TRUE > ${Log}/chr${A}.${Sample}.LDpruned.total.log 2>&1 
done

