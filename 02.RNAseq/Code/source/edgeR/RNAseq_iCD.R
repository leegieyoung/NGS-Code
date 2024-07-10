suppressMessages({
library(edgeR)
library(limma)
library(RColorBrewer)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(repr)
library(statmod)
library(GO.db)
library(ComplexHeatmap)
library(mixOmics)
})

pre_dir <- "/data/keeyoung/gy_RNA/06.output/"

#iCD <- read.table("/data/keeyoung/scRNA/iCD/output/sctype/cibersortx/iCD_count.txt", sep='\t', head=T)
iCD <- read.table("/data/keeyoung/scRNA/iCD/output/sctype/All_Celltype_counts_220923_addBcell.txt", sep='\t') #220923
d <- grep("merge", colnames(iCD))
iCD <- iCD[,-d]
#source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_cpm.R"))
source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_normalized.R"))
RNAseq("iCD",iCD)
#source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_symbol_rpkm.R"))
#RNAseq_rpkm("iCD",iCD)
