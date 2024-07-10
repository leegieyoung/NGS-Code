#!/bin/Rscript
library(ggplot2)
SAMPLE <- Sys.getenv("SAMPLE")
RESULT_DIR <- Sys.getenv("RESULT_DIR")
#setwd(paste0("/ichrogene/project/temp/gylee/0.GWAS/2.plink_result/",SAMPLE))
setwd(RESULT_DIR)

PCA <- read.table("covariate_forplink_PCA.txt", head=T, comment.char="")
Pheno <- read.table(paste0(SAMPLE,"_pheno.txt"), head=T)
Pheno <- Pheno[match(PCA$X.IID, Pheno$IID),]

png(filename=paste0("pca.png"), width=20, height=15, units="cm", res=200)
p<- ggplot(PCA, aes(x=PC1, y=PC2)) +
	geom_point(size=3)

print(p)
dev.off()

PCA <- read.table("covariate_forplink_PCA_outlier.txt", head=T, comment.char="")
Pheno <- read.table(paste0(SAMPLE,"_pheno.txt"), head=T)
Pheno <- Pheno[which(Pheno$IID %in% PCA$X.IID),]

png(filename=paste0("pca_rmoutlier.png"), width=20, height=15, units="cm", res=200)
p<- ggplot(PCA, aes(x=PC1, y=PC2))+
        geom_point(size=3) 
print(p)
dev.off()

