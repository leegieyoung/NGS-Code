#!/bin/Rscript
library(ggplot2);
SAMPLE <- Sys.getenv("SAMPLE")
RESULT_DIR <- Sys.getenv("RESULT_DIR")

print({"---"})
print({"PCA"})
print({"---"})

AGE <- read.table("/mnt/nas/gylee/0.GWAS/1.Input/KNIH_72295ea_info0.8/1.QC/KNIH_72295ea_info0.8.for_age.txt", header=T, comment.char = "")
SEX <- read.table("/mnt/nas/gylee/0.GWAS/1.Input/KNIH_72295ea_info0.8/1.QC/KNIH_72295ea_info0.8.for_sex.txt", header=T, comment.char = "")
head(SEX)
PC20 <- read.table(paste0(RESULT_DIR,SAMPLE,"_PCA.eigenvec"), header=T, comment.char = "")
head(PC20)
print({dim(PC20)})
print({dim(AGE)})
print({dim(SEX)})

mAGE <- match(PC20$X.IID, AGE$IID)
mSEX <- match(PC20$X.IID, SEX$X.IID)
AGE <- AGE$Age[mAGE]
SEX <- SEX$SEX[mSEX]
PHENO <- read.table(paste0(RESULT_DIR,SAMPLE,"_pheno.txt"), header=T, comment.char = "")
wPHENO <- which(PHENO$IID %in% PC20$X.IID)
PHENO <- PHENO$Pheno[wPHENO]
PHENO <- gsub("1","0",PHENO)
PHENO <- gsub("2","1",PHENO)
print({dim(AGE)})
print({dim(SEX)})

COVAR <- cbind(PC20, AGE,SEX, PHENO)
print({head(COVAR)})
colnames(COVAR) <- gsub("X.IID","IID",colnames(COVAR))
write.table(COVAR, paste0(RESULT_DIR,"covariate.txt"), row.names=F, col.names=T, sep='\t', quote=F)

