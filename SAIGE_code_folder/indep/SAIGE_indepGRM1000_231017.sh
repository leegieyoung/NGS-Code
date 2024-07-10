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
Dir="/ichrogene/project/temp/gylee/0.GWAS"
impute_dir="${Dir}/1.Input/${Sample}/0.Impute"
QC_dir="${Dir}/1.Input/${Sample}/1.QC"
MaMi_dir="${Dir}/1.Input/${Sample}/2.MaMi"
Case_pheno="${ana_dir}/temp/11.Case_pheno.txt"
Control_pheno="${ana_dir}/temp/11.Control_pheno.txt"
inversion="${Dir}/REFERENCE/inversion.txt"
Rcode="/ichrogene/project/temp/gylee/0.GWAS/Code/R"
p2_Dir="/ichrogene/project/temp/gylee/Singularity/plink2"
Output_dir="${Dir}/2.plink_result/${Sample}"
#=================================================================
ana_dir="${Dir}/2.plink_result/NoImputed_${Sample}"
SAGE_CD="/ichrogene/project/temp/gylee/0.GWAS/Code/SAIGE_code_folder"


echo "=========================="
echo "Start Step0.createSparseGRM.sh"
echo "=========================="
mkdir -p ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LOG
mkdir -p ${Dir}/2.SAIGE_result/${Sample}/result
Log="${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LOG"


#Made by prune data
#Rscript ${SAGE_CD}/createSparseGRM.R \
#     --plinkFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC/chr${Chr}.LD_pruned_QC_${Sample} \
#     --nThreads=32  \
#     --outputPrefix=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/sparseGRM \
#     --numRandomMarkerforSparseKin=1000 \
#     --relatednessCutoff=0.125 > ${Log}/Stop0.chr${Chr}.log 2>&1
#sleep 3

echo "=========================="
echo "Start step1_fitNULLGLMM.R"
echo "=========================="
mkdir ${Log}/step1_fitNULLGLMM
Rscript ${SAGE_CD}/step1_fitNULLGLMM.R     \
       --sparseGRMFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${Chr}.sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx \
        --sparseGRMSampleIDFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${Chr}.sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
        --plinkFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC/chr${Chr}.LD_pruned_QC_${Sample}  \
        --phenoFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/${Sample}.covariate.txt \
       --useSparseGRMtoFitNULL=TRUE \
       --skipVarianceRatioEstimation=FALSE \
        --phenoCol=PHENO \
       --nThreads=${Thread} \
        --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,AGE,SEX \
	--useSparseGRMforVarRatio=TRUE \
        --qCovarColList=SEX \
        --sampleIDColinphenoFile=IID \
        --traitType=binary \
        --outputPrefix=/ichrogene/project/temp/gylee/0.GWAS/2.SAIGE_result/${Sample}/sparseGRM/fitNGLMM.chr${Chr}.${Sample}.LDpruned.total \
        --IsOverwriteVarianceRatioFile=TRUE > ${Log}/step1_fitNULLGLMM/step1.${Chr}.log 

sleep 3

echo "=========================="
echo "Start step2_SPAtests.R"
echo "=========================="
mkdir ${Log}/step2_SPAtests
Rscript ${SAGE_CD}/step2_SPAtests.R \
        --bedFile=${Output_dir}/chr${Chr}.imputation_g_m_maf_hwe_bfile.bed \
        --bimFile=${Output_dir}/chr${Chr}.imputation_g_m_maf_hwe_bfile.bim \
        --famFile=${Output_dir}/chr${Chr}.imputation_g_m_maf_hwe_bfile.fam \
        --sparseGRMFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${Chr}.sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx \
        --sparseGRMSampleIDFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${Chr}.sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
        --AlleleOrder=alt-first \
        --chrom=${Chr} \
        --GMMATmodelFile=/ichrogene/project/temp/gylee/0.GWAS/2.SAIGE_result/${Sample}/sparseGRM/fitNGLMM.chr${Chr}.${Sample}.LDpruned.total.rda \
        --varianceRatioFile=/ichrogene/project/temp/gylee/0.GWAS/2.SAIGE_result/${Sample}/sparseGRM/fitNGLMM.chr${Chr}.${Sample}.LDpruned.total.varianceRatio.txt \
        --LOCO=FALSE \
        --is_Firth_beta=TRUE \
	--is_imputed_data=TRUE \
        --pCutoffforFirth=0.125 \
        --is_output_moreDetails=TRUE \
        --SAIGEOutputFile=${Dir}/2.SAIGE_result/${Sample}/result/chr${Chr}.${Sample}.LDpruned.total.SAIGE.test.markers.txt \
        --is_fastTest=TRUE \
        --is_output_moreDetails=TRUE > ${Log}/step2_SPAtests/step2.${Chr}.log 
