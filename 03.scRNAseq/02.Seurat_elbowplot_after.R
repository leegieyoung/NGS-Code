#https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
library(Seurat)
library(patchwork)
library(metap)
library(limma)
library(multtest)
library(ggplot2)
library(cowplot)
library(future.apply)
library(dplyr)

pre_dir <- '/data/keeyoung/scRNA/iCD/'
iCD.combined <- readRDS(paste0(pre_dir, "output/1.QC/","raw.iCD.combined" , "_rds"))

iCD.combined <- FindNeighbors(iCD.combined, reduction = "pca", dims = 1:13)
iCD.combined <- FindClusters(iCD.combined, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0) )
iCD.combined <- RunUMAP(object = iCD.combined, dims = 1:13)

cell.ids <- c("C69","C128","C138","C181","C187","C158", "C190", "C193", "C196", "C68", "C129", "C135", "C180", "C186", "C159", "C189", "C192", "C195")
Sample_list <- c("GSM3972009_69","GSM3972013_128","GSM3972016_138","GSM3972020_181","GSM3972022_187","GSM3972017_158","GSM3972024_190","GSM3972026_193","GSM3972028_196", "GSM3972010_68", "GSM3972014_129","GSM3972015_135", "GSM3972019_180", "GSM3972021_186", "GSM3972018_159", "GSM3972023_189", "GSM3972025_192", "GSM3972027_195")

iCD.combined$patients <- '0'
for (i in 1:length(cell.ids)){
print(i)
iCD.combined$patients[grep(paste0(cell.ids[i],"_*"),names(iCD.combined$orig.ident))] <- Sample_list[i]
}

iCD.combined$GIMATs <- '0'
GIMAT <- c("GIMAT","GIMAT","GIMAT","GIMAT","GIMAT","NoGIMAT","NoGIMAT","NoGIMAT","NoGIMAT","Uninf","Uninf","Uninf","Uninf","Uninf","Uninf","Uninf","Uninf","Uninf")
for (i in 1:length(Sample_list)){
    iCD.combined$GIMATs[grep(Sample_list[i], iCD.combined$patients)] <- GIMAT[i]
}

saveRDS(iCD.combined, file = paste0(pre_dir, "output/1.QC/","iCD.combined" , "_rds"))
