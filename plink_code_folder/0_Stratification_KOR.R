#!/bin/Rscript
library(MatchIt)

QC_Dir="/mnt/nas/gylee/0.GWAS/1.Input/KNIH_72295ea_info0.8/1.QC/"

INPUT_Dir=paste0(QC_Dir)
#
test <- read.csv("/mnt/nas/gylee/0.GWAS/1.Input/KNIH_72295ea_info0.8/1.QC/T2D_GRS_age_gender_forR.csv")

test$PHENO <- gsub("Case",1,test$PHENO)
test$PHENO <- gsub("Control",0,test$PHENO)    
test$Gender <- gsub("Male",1,test$Gender)
test$Gender <- gsub("Female",2,test$Gender)

test$PHENO <- as.numeric(test$PHENO)
test$Gender <- as.numeric(test$Gender)

ID <- (sort(test$ID))
m <- match(ID, test$ID)
test <- test[m,]

#Only Gender
test_m <- matchit(PHENO ~ Gender + Age, ratio=2, method="nearest", data=test)
print({test_m})
matched_test <- match.data(test_m)
#print({matched_test})
Case_test <- matched_test[which(matched_test$PHENO==1),]
print({summary(Case_test)})
Control_test <- matched_test[which(matched_test$PHENO==0),]
print({summary(Control_test)})

matched_test <- matched_test[,c("ID","PHENO")]
colnames(matched_test) <- c("IID","PHENO")
matched_test$PHENO <- gsub(1,2,matched_test$PHENO)
matched_test$PHENO <- gsub(0,1,matched_test$PHENO)
IID <- matched_test$IID
write.table(IID ,paste0(INPUT_Dir,"KNIH_72295ea_info0.8.IID"), col.names=T, row.names=F, quote=F)
write.table(matched_test,paste0(INPUT_Dir,"KNIH_72295ea_info0.8","_pheno.txt"), col.names=T, row.names=F, quote=F, sep='\t')
