#made by Giyoung Lee, If you want to help me,contact gy48085@gmail.com
library(Seurat)
library(patchwork)
library(metap)
library(limma)
library(multtest)
source_dir="/data/keeyoung/scRNA/iCD/Code/source/"
source(paste0(source_dir,"Read10X.R"))

pre_dir <- '/data/keeyoung/scRNA/iCD/'
#infla <- c("GSM3972009_69","GSM3972013_128","GSM3972016_138","GSM3972020_181","GSM3972022_187")
#uninfla <- c("GSM3972017_158","GSM3972026_193","GSM3972028_196")
infla <- c("GSM3972024_190","GSM3972030_209")
uninfla <- c("GSM3972010_68","GSM3972014_129","GSM3972015_135","GSM3972018_159","GSM3972019_180","GSM3972021_186","GSM3972023_189","GSM3972025_192","GSM3972027_195","GSM3972029_208")
for (i in 1:length(infla)) {
   dir <- paste0(pre_dir, infla[i])
   read10x(infla[i], dir)
}
   
for (i in 1:length(uninfla)) {
        dir <- paste0(pre_dir, uninfla[i])
        read10x(uninfla[i], dir)
}


