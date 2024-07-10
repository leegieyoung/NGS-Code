#!/bin/Rscript
library(ggplot2);
#PCscore <- read.table("/mnt/nas/gylee/0.GWAS/2.plink_result/NoImputed_Ukb_obesity/PCA/PCA.txt", sep=' ',head=T)

PCA <- read.table("/mnt/nas/gylee/0.GWAS/2.plink_result/Imputed_Ukb_82242_info0.8_stroke/Imputed_Ukb_82242_info0.8_stroke_PCA.eigenvec",head=T, comment.char="")
pheno <- read.table("/mnt/nas/gylee/0.GWAS/2.plink_result/Imputed_Ukb_82242_info0.8_stroke/Imputed_Ukb_82242_info0.8_stroke_pheno.txt", head=T)
w <- which(pheno$IID %in% PCA$X.IID)
PCA$pheno <- pheno$Pheno[w]
PCA$pheno <- gsub(2,"Case",PCA$pheno)
PCA$pheno <- gsub(1,"Control",PCA$pheno)


png(filename="/mnt/nas/gylee/0.GWAS/2.plink_result/Imputed_Ukb_82242_info0.8_stroke/PCA.png", width=20, height=15, units="cm", res=200)
p <- ggplot(PCA, aes(x=PC1, y=PC2, color=pheno))
p <- p + geom_point()
#p < p + scale_color_manual(values=c('#56B4E9','#E69F00'))
p <- p + coord_cartesian(xlim = c(min(PCA$PC1)-0.0001, max(PCA$PC1)+0.0001))
p <- p + coord_cartesian(ylim = c(min(PCA$PC2)-0.0001, max(PCA$PC2)+0.0001))
print({p})
dev.off()
