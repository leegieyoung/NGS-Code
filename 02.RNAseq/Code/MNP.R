source("/data/keeyoung/gy_RNA/Code/RNAseq.R")
MNP <- read.table("/data/keeyoung/scRNA/iCD/output/convertFeatureCounts/MNP_9patients.txt", sep=' ', head=T)

RNAseq("MNP",MNP)
