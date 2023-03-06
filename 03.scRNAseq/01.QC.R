#made by Giyoung Lee, If you want to help me,contact gy48085@gmail.com
library(multtest)
library(ggplot2)
library(cowplot)
library(future.apply)
library(dplyr)

source_dir="/data/keeyoung/scRNA/iCD/Code/source/"
source(paste0(source_dir,"source_Read10X.R"))
source(paste0(source_dir,"source_QC.R"))
pre_dir <- '/data/keeyoung/scRNA/iCD/'

Infla <- c(
"GSM3972009_69",
"GSM3972013_128",
"GSM3972016_138",
"GSM3972020_181",
"GSM3972022_187",
"GSM3972017_158",
"GSM3972024_190",
"GSM3972026_193",
"GSM3972028_196")
#--------------QC
for (i in 1:length(Infla)) {
    QC(pre_dir,Infla[i],"Infla")
}
Uninf <- c(
"GSM3972010_68",
"GSM3972014_129",
"GSM3972015_135",
"GSM3972019_180",
"GSM3972021_186",
"GSM3972018_159",
"GSM3972023_189",
"GSM3972025_192",
"GSM3972027_195")

for (i in 1:length(Uninf)) {
    QC(pre_dir,Uninf[i],"Uninf")
}

