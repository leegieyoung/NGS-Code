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

mkdir ${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LOG
Log="${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LOG/"
#Made by prune data
#Rscript ${SAGE_CD}/createSparseGRM.R       \
#     --plinkFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/LD_pruned_QC_${Sample} \
#     --nThreads=32  \
#     --outputPrefix=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/sparseGRM \
#     --numRandomMarkerforSparseKin=2000 \
#     --relatednessCutoff=0.05

#for A in $(seq 1 22)
#do
#Rscript ${SAGE_CD}/step1_fitNULLGLMM.R     \
#	--sparseGRMFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/sparseGRM_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx \
#        --sparseGRMSampleIDFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/sparseGRM_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
#        --plinkFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${A}.LD_pruned_QC_${Sample}  \
#        --phenoFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/${Sample}.covariate.txt \
#	--useSparseGRMtoFitNULL=TRUE \
#	--skipVarianceRatioEstimation=FALSE \
#        --phenoCol=Pheno \
#	--nThreads=${Thread} \
#        --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,Age,Gender \
#        --qCovarColList=Gender \
#        --sampleIDColinphenoFile=IID \
#        --traitType=binary \
#        --outputPrefix=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/${Sample}/sparseGRM/fitNGLMM.chr${A}.${Sample} \
#        --IsOverwriteVarianceRatioFile=TRUE > ${Log}/chr${A}.${Sample}.log 2>&1
#done
#

#SV-Step2
#for A in $(seq 1 22)
#do
#nohup Rscript ${SAGE_CD}/step2_SPAtests.R \
#        --bedFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${A}.LD_pruned_QC_${Sample}.bed \
#        --bimFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${A}.LD_pruned_QC_${Sample}.bim \
#        --famFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${A}.LD_pruned_QC_${Sample}.fam \
#        --sparseGRMFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/sparseGRM_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx \
#        --sparseGRMSampleIDFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/sparseGRM_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
#        --AlleleOrder=alt-first \
#	--chrom=${A} \
#        --GMMATmodelFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/${Sample}/sparseGRM/fitNGLMM.chr${A}.${Sample}.rda \
#        --varianceRatioFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/${Sample}/sparseGRM/fitNGLMM.chr${A}.${Sample}.varianceRatio.txt \
#        --LOCO=FALSE \
#        --is_Firth_beta=TRUE \
#        --pCutoffforFirth=0.05 \
#        --is_output_moreDetails=TRUE \
#        --SAIGEOutputFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${A}.${Sample}.SAIGE.test.markers.txt \
#        --is_fastTest=TRUE \
#        --is_output_moreDetails=TRUE > ${Log}/chr${A}.${Sample}.SAIGE.test.markers.log 2>&1 &
#done

for A in $(seq 1 22)
do
nohup Rscript ${SAGE_CD}/step1_fitNULLGLMM.R     \
       --sparseGRMFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/sparseGRM_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx \
        --sparseGRMSampleIDFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/sparseGRM_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
        --plinkFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${A}.LD_pruned_QC_${Sample}  \
        --phenoFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/${Sample}.covariate.txt \
       --useSparseGRMtoFitNULL=TRUE \
       --skipVarianceRatioEstimation=FALSE \
        --phenoCol=Pheno \
       --nThreads=${Thread} \
        --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,Age,Gender \
        --qCovarColList=Gender \
        --sampleIDColinphenoFile=IID \
        --traitType=binary \
        --outputPrefix=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/${Sample}/sparseGRM/fitNGLMM.chr${A}.${Sample}.LDpruned.total \
        --IsOverwriteVarianceRatioFile=TRUE > ${Log}/chr${A}.${Sample}.LDpruned.total.log 2>&1 &
done

 
##SV-Step2
#for A in $(seq 1 22)
#do
#nohup Rscript ${SAGE_CD}/step2_SPAtests.R \
#        --bedFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${A}.LD_pruned_QC_${Sample}.bed \
#        --bimFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${A}.LD_pruned_QC_${Sample}.bim \
#        --famFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${A}.LD_pruned_QC_${Sample}.fam \
#        --sparseGRMFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/sparseGRM_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx \
#        --sparseGRMSampleIDFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/sparseGRM_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
#        --AlleleOrder=alt-first \
#       --chrom=${A} \
#        --GMMATmodelFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/${Sample}/sparseGRM/fitNGLMM.chr${A}.${Sample}.LDpruned.total.rda \
#        --varianceRatioFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/${Sample}/sparseGRM/fitNGLMM.chr${A}.${Sample}.LDpruned.total.varianceRatio.txt \
#        --LOCO=FALSE \
#        --is_Firth_beta=TRUE \
#        --pCutoffforFirth=0.05 \
#        --is_output_moreDetails=TRUE \
#        --SAIGEOutputFile=${Dir}/2.SAIGE_result/${Sample}/sparseGRM/chr${A}.${Sample}.LDpruned.total.SAIGE.test.markers.txt \
#        --is_fastTest=TRUE \
#        --is_output_moreDetails=TRUE > ${Log}/chr${A}.${Sample}.LDpruned.total.SAIGE.test.markers.log 2>&1 &
#done
#
#
