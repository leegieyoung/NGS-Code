source("/data/keeyoung/gy_RNA/Code/RNAseq.R")
Tcells <- read.table("/data/keeyoung/scRNA/iCD/output/convertFeatureCounts/Tcells_9patients.txt", sep=' ', head=T)

RNAseq("Tcells",Tcells)
