#!/bin/sh
if [ $# -ne 3 ];then
        echo "Please enter Sample_Name, thread(s)"
               exit
fi
#Sample = AD_14 등으로 기록되어있어야 함
#beagle 진행시Same position 제거를 안했을 경우 사용하는 코드
Sample=$1
Thread=$2
Chr=$3
Dir="/mnt/nas/gylee/0.GWAS"
impute_dir="${Dir}/1.Input/${Sample}/0.Impute"
QC_dir="${Dir}/1.Input/${Sample}/1.QC"
MaMi_dir="${Dir}/1.Input/${Sample}/2.MaMi"
Case_pheno="${ana_dir}/temp/11.Case_pheno.txt"
Control_pheno="${ana_dir}/temp/11.Control_pheno.txt"
inversion="${Dir}/REFERENCE/inversion.txt"
Rcode="/mnt/nas/gylee/0.GWAS/Code/R"
p2_Dir="/mnt/nas/gylee/Singurality/plink2"
Output_dir="${Dir}/2.plink_result/${Sample}"
#=================================================================
ana_dir="${Dir}/2.plink_result/NoImputed_${Sample}"
SAGE_CD="/mnt/nas/gylee/0.GWAS/Code/SAIGE_code_folder"


echo "=========================="
echo "Start Step0.createSparseGRM.sh"
echo "=========================="
mkdir -p ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LOG
mkdir -p ${Dir}/2.SAIGE_result/${Sample}/result
Log="${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LOG"


#Made by prune data
Rscript ${SAGE_CD}/createSparseGRM.R \
     --plinkFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC/chr${Chr}.LD_pruned_QC_${Sample} \
     --nThreads=32  \
     --outputPrefix=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${Chr}.sparseGRM \
     --numRandomMarkerforSparseKin=1000 \
     --relatednessCutoff=0.125 > ${Log}/Stop0.chr${Chr}.log 2>&1
sleep 3

echo "=========================="
echo "Start step1_fitNULLGLMM.R"
echo "=========================="
Rscript ${SAGE_CD}/step1_fitNULLGLMM.R     \
       --sparseGRMFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${Chr}.sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx \
        --sparseGRMSampleIDFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${Chr}.sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
        --plinkFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC/chr${Chr}.LD_pruned_QC_${Sample}  \
        --phenoFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/${Sample}.covariate.txt \
       --useSparseGRMtoFitNULL=TRUE \
       --skipVarianceRatioEstimation=FALSE \
        --phenoCol=PHENO \
	--LOCO=FALSE \
       --nThreads=${Thread} \
        --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,AGE,SEX \
        --qCovarColList=SEX \
        --sampleIDColinphenoFile=IID \
        --traitType=binary \
        --outputPrefix=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/${Sample}/sparseGRM/fitNGLMM.chr${Chr}.${Sample}.LDpruned.total_1000 \
        --IsOverwriteVarianceRatioFile=TRUE > ${Log}/chr${Chr}.${Sample}.LDpruned.total_1000.log 2>&1

sleep 3

echo "=========================="
echo "Start step2_SPAtests.R"
echo "=========================="
Rscript ${SAGE_CD}/step2_SPAtests.R \
        --bedFile=${Output_dir}/chr${Chr}.imputation_g_m_maf_hwe_bfile.bed \
        --bimFile=${Output_dir}/chr${Chr}.imputation_g_m_maf_hwe_bfile.bim \
        --famFile=${Output_dir}/chr${Chr}.imputation_g_m_maf_hwe_bfile.fam \
        --sparseGRMFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${Chr}.sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx \
        --sparseGRMSampleIDFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${Chr}.sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
        --AlleleOrder=alt-first \
        --chrom=${Chr} \
        --GMMATmodelFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/${Sample}/sparseGRM/fitNGLMM.chr${Chr}.${Sample}.LDpruned.total_1000.rda \
        --varianceRatioFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/${Sample}/sparseGRM/fitNGLMM.chr${Chr}.${Sample}.LDpruned.total_1000.varianceRatio.txt \
        --LOCO=FALSE \
        --is_Firth_beta=TRUE \
        --pCutoffforFirth=0.125 \
        --is_output_moreDetails=TRUE \
        --SAIGEOutputFile=${Dir}/2.SAIGE_result/${Sample}/result/chr${Chr}.${Sample}.LDpruned.total.SAIGE.test.markers_1000.txt \
        --is_fastTest=TRUE \
        --is_output_moreDetails=TRUE > ${Log}/chr${Chr}.${Sample}.LDpruned.total.SAIGE.test.markers_1000.log 3>&1
