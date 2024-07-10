source("/data/keeyoung/gy_RNA/Code/RNAseq.R")
Plasma <- read.table("/data/keeyoung/scRNA/iCD/output/convertFeatureCounts/Plasma_9patients.txt", sep=' ', head=T)

RNAseq("Plasma",Plasma)
