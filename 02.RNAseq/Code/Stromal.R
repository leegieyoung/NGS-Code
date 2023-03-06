source("/data/keeyoung/gy_RNA/Code/RNAseq.R")
Stromal <- read.table("/data/keeyoung/scRNA/iCD/output/convertFeatureCounts/Stromal_9patients.txt", sep=' ', head=T)

RNAseq("Stromal",Stromal)
