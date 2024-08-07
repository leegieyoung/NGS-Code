#!/bin/Rscript
library(ggplot2);
SAMPLE <- Sys.getenv("SAMPLE")
RESULT_DIR <- Sys.getenv("RESULT_DIR")

print({"---"})
print({"PCA"})
print({"---"})

AGE <- read.table("/mnt/nas/gylee/0.GWAS/1.Input/Ukb_obesity/1.raw/subset_of_ukb45411_for_age.txt", header=T, comment.char = "")
SEX <- read.table("/mnt/nas/gylee/0.GWAS/1.Input/Ukb_obesity/1.raw/subset_of_ukb45411_for_sex.txt", header=T, comment.char = "")
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
print({dim(AGE)})
print({dim(SEX)})

COVAR <- cbind(PC20, AGE,SEX)
colnames(COVAR) <- gsub("X.IID","#IID",colnames(COVAR))
write.table(COVAR, paste0(RESULT_DIR,"covariate_forplink.txt"), row.names=F, col.names=T, sep='\t', quote=F)

