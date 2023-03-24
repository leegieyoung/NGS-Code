#source("/data/keeyoung/gy_RNA/Code/RNAseq.R")
source("RNAseq_symbol_rpkm.R")
#iCD <- read.table("/data/keeyoung/scRNA/iCD/output/convertFeatureCounts/iCD_9patients.txt", sep=' ', head=T)
iCD <- read.table("/data/keeyoung/scRNA/iCD/output/sctype/cibersortx/iCD_count.txt", sep='\t')
RNAseq("iCD",iCD)

