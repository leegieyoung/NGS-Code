Rscript createSparseGRM.R       \
     --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
     --nThreads=4  \
     --outputPrefix=./output/sparseGRM       \
     --numRandomMarkerforSparseKin=2000      \
     --relatednessCutoff=0.125

