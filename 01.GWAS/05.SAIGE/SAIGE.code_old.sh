p2_Dir="/mnt/nas/gylee/Singurality/plink2"
SAGE_CD="/mnt/nas/gylee/0.GWAS/Code/SAIGE_code_folder"


#Made by prune data
Rscript ${SAGE_CD}/createSparseGRM.R       \
     --plinkFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity/sparseGRM/Cli_LD_pruned_QC_Ukb_obesity \
     --nThreads=31  \
     --outputPrefix=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity/sparseGRM/sparseGRM \
     --numRandomMarkerforSparseKin=1000 \
     --relatednessCutoff=0.05

#for A in $(seq 1 22)
#do
#Rscript ${SAGE_CD}/step1_fitNULLGLMM.R     \
#	--sparseGRMFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity/sparseGRM/sparseGRM_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx \
#        --sparseGRMSampleIDFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity/sparseGRM/sparseGRM_relatednessCutoff_0.05_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
#	--useSparseGRMtoFitNULL=TRUE \
#        --plinkFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity/sparseGRM/chr${A}.1000_for_vr_Ukb_obesity  \
#        --phenoFile=/mnt/nas/gylee/0.GWAS/1.Input/Ukb_obesity/3.covaraite/covariate.txt \
#        --phenoCol=Pheno \
#        --covarColList=PC1,PC2,PC3,PC4,Age,Gender \
#        --qCovarColList=Gender \
#        --sampleIDColinphenoFile=IID \
#        --traitType=binary \
#	--isCateVarianceRatio=TRUE \
#        --outputPrefix=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/chr${A}.Ukb_obesity \
#        --invNormalize=FALSE \
#	--memoryChunk=250 \
#        --IsOverwriteVarianceRatioFile=TRUE
#done

##Fail..
##Error in Get_Variance_Ratio(varianceRatioFile, cateVarRatioMinMACVecExclude,  : 
##  sparse GRM is specified but the variance ratio for sparse GRM was not estimatedin Step 1.
#Rscript ${SAGE_CD}/step1_fitNULLGLMM.R     \
#        --sparseGRMFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity/sparseGRM/sparseGRM_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
#        --sparseGRMSampleIDFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity/sparseGRM/sparseGRM_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
#        --plinkFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity/sparseGRM/1000_for_vr_Ukb_obesity  \
#        --phenoFile=/mnt/nas/gylee/0.GWAS/1.Input/Ukb_obesity/3.covaraite/covariate.txt \
#        --qCovarColList=Gender \
#        --phenoCol=Pheno \
#        --covarColList=PC1,PC2,PC3,PC4,Age,Gender \
#        --sampleIDColinphenoFile=IID \
#        --traitType=binary \
#        --outputPrefix=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity \
#        --nThreads=31 \
#        --invNormalize=FALSE \
#        --IsOverwriteVarianceRatioFile=TRUE

##Fail
#Rscript ${SAGE_CD}/step1_fitNULLGLMM.R     \
#	--sparseGRMFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity/sparseGRM/sparseGRM_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
#	--sparseGRMSampleIDFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity/sparseGRM/sparseGRM_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
#        --plinkFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity/sparseGRM/1000_for_vr_Ukb_obesity  \
#	--useSparseGRMtoFitNULL=TRUE \
#        --phenoFile=/mnt/nas/gylee/0.GWAS/1.Input/Ukb_obesity/3.covaraite/covariate.txt \
#	--qCovarColList=Gender \
#	--skipVarianceRatioEstimation=FALSE \
#	--phenoCol=Pheno \
#        --covarColList=PC1,PC2,PC3,PC4,Age,Gender \
#        --sampleIDColinphenoFile=IID \
#        --traitType=binary \
#        --outputPrefix=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity \
#        --nThreads=31 \
#        --IsOverwriteVarianceRatioFile=TRUE

##Based Single_Step2
#nohup Rscript ${SAGE_CD}/step2_SPAtests.R \
#	--bedFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity/sparseGRM/1000_for_vr_Ukb_obesity.bed \
#	--bimFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity/sparseGRM/1000_for_vr_Ukb_obesity.bim \
#	--famFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity/sparseGRM/1000_for_vr_Ukb_obesity.fam \
#        --AlleleOrder=alt-first \
#        --SAIGEOutputFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity/Ukb_obesity_SAIGE.markers.txt \
#	--chrom=1 \
#	--minMAF=0 \
#	--minMAC=20 \
#        --GMMATmodelFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity.rda \
#        --varianceRatioFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity.varianceRatio.txt \
#        --LOCO=FALSE \
#        --is_Firth_beta=TRUE    \
#        --pCutoffforFirth=0.1 \
#        --is_output_moreDetails=TRUE \
#        --sparseGRMFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity/sparseGRM/sparseGRM_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
#	--sparseGRMSampleIDFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity/sparseGRM/sparseGRM_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
#	--is_fastTest=TRUE > /mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity/LOG/Ukb_obesity_SAIGE.markers.log 2>&1 &

#High qual variants.. but fall, uk exome
#nohup Rscript ${SAGE_CD}/step2_SPAtests.R \
#        --bedFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity/sparseGRM/1000_for_vr_Ukb_obesity.bed \
#        --bimFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity/sparseGRM/1000_for_vr_Ukb_obesity.bim \
#        --famFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity/sparseGRM/1000_for_vr_Ukb_obesity.fam \
#        --sparseGRMFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity/sparseGRM/sparseGRM_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
#        --sparseGRMSampleIDFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity/sparseGRM/sparseGRM_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
#        --AlleleOrder=alt-first \
#        --GMMATmodelFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity.rda \
#        --varianceRatioFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity.varianceRatio.txt \
#        --LOCO=FALSE \
#        --is_Firth_beta=TRUE \
#        --pCutoffforFirth=0.05 \
#        --is_output_moreDetails=TRUE \
#        --SAIGEOutputFile=/mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity/Ukb_obesity_SAIGE.markers.txt \
#        --is_fastTest=TRUE \
#        --is_output_moreDetails=TRUE > /mnt/nas/gylee/0.GWAS/2.SAIGE_result/Ukb_obesity/LOG/Ukb_obesity_SAIGE.markers.log 2>&1 &
