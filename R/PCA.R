#!/bin/Rscript
library(ggplot2);
SAMPLE <- Sys.getenv("SAMPLE")
QC_DIR <- Sys.getenv("QC_DIR")

print({"---"})
print({"PCA"})
print({"---"})

PCscore <- read.table(paste0(QC_DIR,"PCA.txt"), header=T, comment.char = "")

png(filename=paste0(QC_DIR,"PCA.png"), width=20, height=15, units="cm", res=200)
p <- ggplot(PCscore, aes(x=PC1, y=PC2, color=Pheno))
p <- p + geom_point()
#p < p + scale_color_manual(values=c('#56B4E9','#E69F00'))
p <- p + coord_cartesian(xlim = c(min(PCscore$PC1)-0.0001, max(PCscore$PC1)+0.0001))
p <- p + coord_cartesian(ylim = c(min(PCscore$PC2)-0.0001, max(PCscore$PC2)+0.0001))
print({p})
dev.off()

AGE <- read.table("/mnt/nas/gylee/0.GWAS/1.Input/Ukb_obesity/1.raw/subset_of_ukb45411_for_age.txt", header=T, comment.char = "")
SEX <- read.table("/mnt/nas/gylee/0.GWAS/1.Input/Ukb_obesity/1.raw/subset_of_ukb45411_for_sex.txt", header=T, comment.char = "")
head(SEX)
PC20 <- read.table(paste0(QC_DIR,SAMPLE,"_g_m_maf_hwe_het_PCA.eigenvec"), header=T, comment.char = "")
mAGE <- match(PC20$IID, AGE$IID)
mSEX <- match(PC20$IID, SEX$IID)
AGE <- AGE$Age[mAGE]
SEX <- SEX$SEX[mSEX]

COVAR <- cbind(PC20, AGE,SEX)
colnames(COVAR) <- gsub("X.FID","#FID",colnames(COVAR))
write.table(COVAR, paste0(QC_DIR,"covariate.txt"), row.names=F, col.names=T, sep='\t', quote=F)

