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
Output_dir="${Dir}/2.plink_result/${Sample}"

#=================================================================
ana_dir="${Dir}/2.plink_result/NoImputed_${Sample}"
SAGE_CD="/mnt/nas/gylee/0.GWAS/Code/SAIGE_code_folder"

echo "=========================="
echo "Start Step1.fitNULLGLMM.sh"
echo "=========================="

export QC_DIR="${QC_dir}/"
export RESULT_DIR="${Output_dir}/"
echo ${ANA_DIR}
export SAMPLE="${Sample}"


mkdir -p ${Dir}/2.SAIGE_result/temp/sparseGRM/LOG
mkdir -p ${Dir}/2.SAIGE_result/temp/result
Log="${Dir}/2.SAIGE_result/temp/sparseGRM/LOG/"

singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/GWAS.sif Rscript /mnt/nas/gylee/0.GWAS/Code/R/covariate_forSAIGE_KOR.R
mkdir -p ${Dir}/2.SAIGE_result/temp/sparseGRM/LOG
cp ${Output_dir}/covariate.txt ${Dir}/2.SAIGE_result/temp/sparseGRM/${Sample}.covariate.txt

A=6

#Step1
singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/SAIGE/SAIGE_docker_R36.sif Rscript ${SAGE_CD}/step1_fitNULLGLMM.R     \
       --sparseGRMFile=${Dir}/2.SAIGE_result/temp/sparseGRM/chr${A}.sparseGRM_relatednessCutoff_0.125_5000_randomMarkersUsed.sparseGRM.mtx \
        --sparseGRMSampleIDFile=${Dir}/2.SAIGE_result/temp/sparseGRM/chr${A}.sparseGRM_relatednessCutoff_0.125_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
        --plinkFile=${Dir}/2.SAIGE_result/temp/sparseGRM/chr${A}.LD_pruned_QC_${Sample}  \
        --phenoFile=${Dir}/2.SAIGE_result/temp/sparseGRM/${Sample}.covariate.txt \
       --useSparseGRMtoFitNULL=TRUE \
       --skipVarianceRatioEstimation=FALSE \
        --phenoCol=PHENO \
       --nThreads=${Thread} \
        --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,AGE,SEX \
        --qCovarColList=SEX \
        --sampleIDColinphenoFile=IID \
        --traitType=binary \
        --outputPrefix=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/temp/sparseGRM/fitNGLMM.chr${A}.${Sample}.LDpruned.total \
        --IsOverwriteVarianceRatioFile=TRUE > ${Log}/chr${A}.${Sample}.LDpruned.total.log 2>&1 

singularity exec --bind /mnt/:/mnt/ /mnt/nas/gylee/Singurality/SAIGE/SAIGE_docker_R36.sif Rscript ${SAGE_CD}/step2_SPAtests.R \
        --bedFile=${Output_dir}/chr${A}.imputation_g_m_maf_hwe_bfile.bed \
        --bimFile=${Output_dir}/chr${A}.imputation_g_m_maf_hwe_bfile.bim \
        --famFile=${Output_dir}/chr${A}.imputation_g_m_maf_hwe_bfile.fam \
        --sparseGRMFile=${Dir}/2.SAIGE_result/temp/sparseGRM/chr${A}.sparseGRM_relatednessCutoff_0.125_5000_randomMarkersUsed.sparseGRM.mtx \
        --sparseGRMSampleIDFile=${Dir}/2.SAIGE_result/temp/sparseGRM/chr${A}.sparseGRM_relatednessCutoff_0.125_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
        --AlleleOrder=alt-first \
       --chrom=${A} \
        --GMMATmodelFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/temp/sparseGRM/fitNGLMM.chr${A}.${Sample}.LDpruned.total.rda \
        --varianceRatioFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/temp/sparseGRM/fitNGLMM.chr${A}.${Sample}.LDpruned.total.varianceRatio.txt \
        --LOCO=FALSE \
        --is_Firth_beta=TRUE \
        --pCutoffforFirth=0.125 \
        --is_output_moreDetails=TRUE \
        --SAIGEOutputFile=${Dir}/2.SAIGE_result/temp/result/chr${A}.${Sample}.LDpruned.total.SAIGE.test.markers.txt \
        --is_fastTest=TRUE \
        --is_output_moreDetails=TRUE > ${Log}/chr${A}.${Sample}.LDpruned.total.SAIGE.test.markers.log 2>&1
