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
AGE <- AGE[mAGE,]
SEX <- SEX[mSEX,]
PHENO <- read.table(paste0(RESULT_DIR,SAMPLE,"_pheno.txt"), header=T, comment.char = "")
mPHENO <- match(PC20$X.IID, PHENO$IID)
PHENO <- PHENO[mPHENO,]
PHENO$Pheno <- gsub("1","0",PHENO$Pheno)
PHENO$Pheno <- gsub("2","1",PHENO$Pheno)

rownames(AGE) <- AGE$IID
AGE <- AGE[,-which(colnames(AGE) %in% "IID")]
rownames(SEX) <- SEX$X.IID
SEX <- SEX[,-which(colnames(SEX) %in% "X.IID")]
rownames(PHENO) <- PHENO$IID
PHENO <- PHENO[,-which(colnames(PHENO) %in% "IID")]
rownames(PC20) <- PC20$X.IID

#COVAR <- cbind(PC20, AGE,SEX, PHENO)
if(all.equal(rownames(AGE), rownames(SEX), rownames(PC20), rownames(PHENO))=="TRUE"){
	COVAR <- cbind(PC20, AGE,SEX, PHENO)
}
colnames(COVAR) <- gsub("X.IID","IID",colnames(COVAR))
write.table(COVAR, paste0(RESULT_DIR,"covariate.txt"), row.names=F, col.names=T, sep='\t', quote=F)

