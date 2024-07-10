#!/bin/Rscript
library(MatchIt)

QC_Dir="/mnt/nas/gylee/0.GWAS/1.Input/"
Sample="Ukb_T2D_GB"
INPUT_Dir=paste0(QC_Dir,Sample,"/1.QC/")
#
test <- read.table(paste0(INPUT_Dir,"sample/","subset_of_ukb45411_for_AgeSex_",Sample,".txt"), head=T)

test$Gender <- gsub("Female",2,test$Gender)
test$Gender <- gsub("Male",1,test$Gender)
test$PHENO <- gsub("Case",1,test$PHENO)
test$PHENO <- gsub("Control",0,test$PHENO)
test$PHENO <- as.numeric(test$PHENO)
test$Gender <- as.numeric(test$Gender)

#Only Gender
test_m <- matchit(PHENO ~ Gender, ratio=1, method="nearest",data=test)
print({test_m})
matched_test <- match.data(test_m)
#print({matched_test})
Case_test <- matched_test[which(matched_test$PHENO==1),]
print({summary(Case_test)})
Control_test <- matched_test[which(matched_test$PHENO==0),]
print({summary(Control_test)})

matched_test <- matched_test[,c("FID","PHENO")]
colnames(matched_test) <- c("IID","PHENO")
matched_test$PHENO <- gsub(1,2,matched_test$PHENO)
matched_test$PHENO <- gsub(0,1,matched_test$PHENO)
IID <- matched_test$IID
write.table(IID ,paste0(INPUT_Dir,Sample,".IID"), col.names=T, row.names=F, quote=F)
write.table(matched_test,paste0(INPUT_Dir,Sample,"_pheno.txt"), col.names=T, row.names=F, quote=F, sep='\t')
