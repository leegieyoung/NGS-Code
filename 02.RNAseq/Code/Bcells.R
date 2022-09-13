source("/data/keeyoung/gy_RNA/Code/RNAseq.R")
Bcells <- read.table("/data/keeyoung/scRNA/iCD/output/convertFeatureCounts/Bcells_9patients.txt", sep=' ', head=T)

RNAseq("Bcells",Bcells)
