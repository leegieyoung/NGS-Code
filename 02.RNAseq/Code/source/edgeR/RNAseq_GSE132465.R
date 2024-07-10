GSE132465_celltype <- read.table("/data/keeyoung/scRNA/CRC/01.raw.data/scRNA/GSE132465/GSE132465_anno_GY.txt", sep='\t')
GSE132465_celltype <- cbind(rownames(GSE132465_celltype), GSE132465_celltype[,5])
GSE132465_celltype <- gsub("-",".",GSE132465_celltype)
iCD <- read.table("/data/keeyoung/scRNA/CRC/01.raw.data/scRNA/GSE132465/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz", sep='\t')
#iCD <- as.matrix(iCD)
colnames(iCD)<- GSE132465_celltype[,2][match(colnames(iCD), GSE132465_celltype[,1])]

pre_dir <- "/data/keeyoung/gy_RNA/06.output/"
source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_normalized.R"))
RNAseq("GSE132465",iCD)

