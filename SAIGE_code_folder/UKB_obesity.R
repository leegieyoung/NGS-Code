Rscript step1_fitNULLGLMM.R     \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr  \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --phenoCol=y_binary \
        --covarColList=x1,x2 \
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=./output/example_binary \
        --nThreads=1	\
        --IsOverwriteVarianceRatioFile=TRUE

Rscript createSparseGRM.R       \
     --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
     --nThreads=4  \
     --outputPrefix=./output/sparseGRM       \
     --numRandomMarkerforSparseKin=2000      \
     --relatednessCutoff=0.125

#Fit the null model and estimate a variance ratio
Rscript step1_fitNULLGLMM.R     \
        --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx   \
        --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt     \
        --useSparseGRMtoFitNULL=TRUE    \
        --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr \
        --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
        --skipVarianceRatioEstimation=FALSE \
        --phenoCol=y_binary \
        --covarColList=x1,x2 \
        --qCovarColList=x2  \
        --sampleIDColinphenoFile=IID \
        --traitType=binary        \
        --outputPrefix=./output/example_binary_sparseGRM_vr \
        --IsOverwriteVarianceRatioFile=TRUE

##Step 2: 
#performing single-variant association tests
Rscript step2_SPAtests.R        \
        --bgenFile=./input/genotype_100markers.bgen    \
        --bgenFileIndex=./input/genotype_100markers.bgen.bgi \
        --sampleFile=./input/samplelist.txt \
        --AlleleOrder=ref-first \
        --SAIGEOutputFile=./output/genotype_100markers_marker_bgen_fullGRMforNull_with_vr.txt    \
        --chrom=1       \
        --minMAF=0 \
        --minMAC=20 \
        --GMMATmodelFile=./output/example_binary.rda \
        --varianceRatioFile=./output/example_binary.varianceRatio.txt \
        --is_Firth_beta=TRUE    \
        --pCutoffforFirth=0.05 \
        --is_output_moreDetails=TRUE    \
        --LOCO=TRUE

#Sparse GRM was used for fitting the model in Step 1. Variance ratio is estimated
Rscript step2_SPAtests.R        \
        --bgenFile=./input/genotype_100markers.bgen    \
        --bgenFileIndex=./input/genotype_100markers.bgen.bgi \
        --sampleFile=./input/samplelist.txt \
        --AlleleOrder=ref-first \
        --SAIGEOutputFile=./output/genotype_100markers_marker_bgen_Firth.txt    \
        --chrom=1       \
        --minMAF=0 \
        --minMAC=20 \
        --GMMATmodelFile=./output/example_binary_sparseGRM_vr.rda \
        --varianceRatioFile=./output/example_binary_sparseGRM_vr.varianceRatio.txt  \
        --LOCO=FALSE \
        --is_Firth_beta=TRUE    \
        --pCutoffforFirth=0.1 \
        --is_output_moreDetails=TRUE \
        --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx   \
        --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
        --is_fastTest=TRUE
