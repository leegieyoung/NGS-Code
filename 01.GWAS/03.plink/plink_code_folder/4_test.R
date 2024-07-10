SAMPLE <- Sys.getenv("SAMPLE")
DIR <- Sys.getenv("DIR")
ANA_DIR <- Sys.getenv("ANA_DIR")

print(paste0("SAMPLE is ",SAMPLE))
print(paste0("DIR is ",DIR))
print(paste0("ANA_DIR is ", ANA_DIR))
agesex <- read.csv("/mnt/nas/gylee/0.GWAS/2.SAIGE_result/subset_of_ukb45411_for_AgeSex.csv")
#PCscore <- read.table("/mnt/nas/gylee/0.GWAS/2.plink_result/NoImputed_Ukb_obesity/PCA/Ukb_obesity_NoNA_PCA.eigenvec", head=T, comment.char = "")
PCscore <- read.table(paste0("/mnt/nas/gylee/0.GWAS/2.plink_result/NoImputed_",SAMPLE,"/PCA/",SAMPLE,"_NoNA_PCA.eigenvec"), head=T, comment.char = "")

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

#fam <- read.table("/mnt/nas/gylee/0.GWAS/1.Input/Ukb_obesity/2.MaMi/QC_Imputed_Ukb_obesity.fam", head=F)
#mf <- match(covariate$IID, fam$V1)
#covariate  <- cbind(covariate, fam$V6[mf])
#covariate <- covariate[c("IID","PC1","PC2","PC3","PC4","Age","Gender","fam$V6[mf]")]
#colnames(covariate) <- c("IID","PC1","PC2","PC3","PC4","Age","Gender","Pheno")
#which(is.na(covariate$Pheno))
#covariate$Pheno <- gsub(2,"Case",covariate$Pheno)
#covariate$Pheno <- gsub(1,"Control",covariate$Pheno)
#covariate$Pheno <- gsub("Case","1",covariate$Pheno)
#covariate$Pheno <- gsub("Control","0",covariate$Pheno)

#write.table(covariate, paste0(DIR,SAMPLE,"covariate.txt"), sep='\t', quote=F, col.names=T, row.names=F)

fam <- read.table(paste0(ANA_DIR,"temp/8.raw_",SAMPLE,"_NoNA.fam"), head=F)
covariate <- covariate[which(covariate$IID %in% fam$V1),]
fam <- fam[which(fam$V1 %in% covariate$IID),]
mf <- match(covariate$IID, fam$V1)

paste0("FamID and CovariateID is all equal? ",all.equal(fam[mf,"V1"], covariate$IID))
covariate  <- cbind(covariate, fam$V6[mf])
covariate <- covariate[c("IID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","Age","Gender","fam$V6[mf]")]
colnames(covariate) <- c("IID","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20","Age","Gender","Pheno")
which(is.na(covariate$Pheno))
covariate$Pheno <- gsub(2,"Case",covariate$Pheno)
covariate$Pheno <- gsub(1,"Control",covariate$Pheno)
covariate$Pheno <- gsub("Case","1",covariate$Pheno)
covariate$Pheno <- gsub("Control","0",covariate$Pheno)

write.table(covariate, paste0(DIR,SAMPLE,".covariate.txt"), sep='\t', quote=F, col.names=T, row.names=F)

