agesex <- read.csv("/mnt/nas/gylee/0.GWAS/2.SAIGE_result/subset_of_ukb45411_for_AgeSex.csv")
PCscore <- read.table("/mnt/nas/gylee/0.GWAS/2.plink_result/NoImputed_Ukb_obesity/PCA/Ukb_obesity_NoNA_PCA.eigenvec", head=T, comment.char = "")

wPC <- which(PCscore$IID %in% agesex$FID)
PCscore <- PCscore[wPC,]
wax <- which(agesex$FID %in% PCscore$IID)
agesex <- agesex[wax,]

m <- match(PCscore$IID, agesex$FID)
covariate <- cbind(PCscore, agesex[m,])
which(is.na(covariate$Gender))
which(is.na(covariate$Age))

covariate$Gender <- gsub("Male",2,covariate$Gender)
covariate$Gender <- gsub("Female",1,covariate$Gender)

fam <- read.table("/mnt/nas/gylee/0.GWAS/1.Input/Ukb_obesity/2.MaMi/QC_Imputed_Ukb_obesity.fam", head=F)
mf <- match(covariate$IID, fam$V1)

covariate  <- cbind(covariate, fam$V6[mf])
covariate <- covariate[c("IID","PC1","PC2","PC3","PC4","Age","Gender","fam$V6[mf]")]
colnames(covariate) <- c("IID","PC1","PC2","PC3","PC4","Age","Gender","Pheno")
which(is.na(covariate$Pheno))
covariate$Pheno <- gsub(2,"Case",covariate$Pheno)
covariate$Pheno <- gsub(1,"Control",covariate$Pheno)
covariate$Pheno <- gsub("Case","1",covariate$Pheno)
covariate$Pheno <- gsub("Control","0",covariate$Pheno)

write.table(covariate, "/mnt/nas/gylee/0.GWAS/1.Input/Ukb_obesity/3.covaraite/covariate.txt", sep='\t', quote=F, col.names=T, row.names=F)
