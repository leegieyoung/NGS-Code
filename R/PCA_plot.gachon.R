#!/bin/Rscript
library(ggplot2)
setwd("/ichrogene/project/temp/gylee/0.GWAS/2.plink_result/gachon")

PCA <- read.table("covariate_forplink_PCA.txt", head=T, comment.char="")
Pheno <- read.table("gachon_pheno.txt", head=T)
Pheno[match(PCA$X.IID, Pheno$IID),]
PCA$Pheno <- Pheno[match(PCA$X.IID, Pheno$IID),"Pheno"]
PCA$Pheno <- gsub(1,"Control",PCA$Pheno)
PCA$Pheno <- gsub(2,"Case",PCA$Pheno)

png(filename=paste0("pca.png"), width=20, height=15, units="cm", res=200)
p<- (ggplot(PCA, aes(x=PC1, y=PC2, colour=factor(Pheno), shape=factor(Pheno))) +
	geom_point(size=3) +
	scale_color_manual(values=c('#56B4E9','#E69F00')) +
        labs(color = "Phenotype", shape = "Phenotype")
)
print(p)
dev.off()

PCA <- read.table("covariate_forplink_PCA_outlier.txt", head=T, comment.char="")
Pheno <- read.table("gachon_pheno.txt", head=T)
Pheno <- Pheno[which(Pheno$IID %in% PCA$IID),]
Pheno[match(PCA$X.IID, Pheno$IID),]
PCA$Pheno <- Pheno[match(PCA$IID, Pheno$IID),"Pheno"]
PCA$Pheno <- gsub(1,"Control",PCA$Pheno)
PCA$Pheno <- gsub(2,"Case",PCA$Pheno)

png(filename=paste0("pca_rmoutlier.png"), width=20, height=15, units="cm", res=200)
p<- (ggplot(PCA, aes(x=PC1, y=PC2, colour=factor(Pheno), shape=factor(Pheno))) +
        geom_point(size=3) +
        scale_color_manual(values=c('#56B4E9','#E69F00')) +
	labs(color = "Phenotype", shape = "Phenotype")
)
print(p)
dev.off()

