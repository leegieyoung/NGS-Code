library(MatchIt)
WC <- read.table("/mnt/nas/gylee/0.GWAS/Ukb/48/48.result", head=T)
WC$result <- WC[,1]
for (A in 2:length(colnames(WC))){
	w <- which(is.na(WC[,A]))
	WC$result[-w] <- WC[-w,A]
}
WC_result <- as.data.frame(cbind(WC$eid,WC$result))
colnames(WC_result) <- c("eid","cm")
WC_result$inch <- WC_result$cm/2.54
WC_result$eid <- as.character(WC_result$eid)
#write.table(WC_result, "/mnt/nas/gylee/0.GWAS/Ukb/21001/WC.result", col.names=T, row.names=F, quote=F)

BMI <- read.table("/mnt/nas/gylee/0.GWAS/Ukb/21001/subset_of_ukb45411_for_BMI.csv", head=T, sep=',', fill=TRUE)
BMI$result <- BMI[,2]
for (A in 3:length(colnames(BMI))){
        w <- which(is.na(BMI[,A]))
        length(w)
        BMI$result[-w] <- BMI[-w,A]
}
BMI_result <- as.data.frame(cbind(BMI$eid,BMI$result))
colnames(BMI_result) <- c("eid","BMI")
BMI_result$eid <- as.character(BMI_result$eid)
#write.table(BMI_result, "/mnt/nas/gylee/0.GWAS/Ukb/21001/BMI.result", col.names=T, row.names=F, quote=F)
w1 <- which(WC_result$eid %in% BMI_result$eid)
all.equal(BMI_result$eid, WC_result$eid[w1])
BMI_result$inch <- WC_result$inch[w1]

AgeSex <- read.table("/mnt/nas/gylee/0.GWAS/Ukb/subset_of_ukb45411_for_AgeSex.txt", head=T)
AgeSex$FID <- as.character(AgeSex$FID)
wAgeSex <- which(AgeSex$FID %in% BMI_result$eid)
AgeSex <- AgeSex[wAgeSex,]
wBMI_result <- which(BMI_result$eid %in% AgeSex$FID)
BMI_result <- BMI_result[wBMI_result,]
all.equal(BMI_result$eid, AgeSex$FID)
BMI_result$Age <- AgeSex$Age
BMI_result$Gender <- AgeSex$Gender

rmlist <- which(is.na(as.numeric(BMI_result$Age)))
BMI_result <- BMI_result[-rmlist,]
BMI_result$Age <- as.numeric(BMI_result$Age)

BMI_result$Gender <- gsub("Male",1,BMI_result$Gender)
BMI_result$Gender <- gsub("Female",2,BMI_result$Gender)
BMI_result$Gender <- as.numeric(BMI_result$Gender)

summary(BMI_result)
if(length(which(BMI_result$BMI >= 100)) !=0){
	BMI_result <- BMI_result[-which(BMI_result$BMI >= 100),]
}
if(length(which(BMI_result$inch >=200)) !=0){
	BMI_result <- BMI_result[-which(BMI_result$inch >= 200),]
}
if(length(which(BMI_result$inch < 15)) !=0){
	BMI_result <- BMI_result[-which(BMI_result$inch < 15),]
}
if(length(which(is.na(BMI_result$BMI))) !=0){
        BMI_result <- BMI_result[-which(is.na(BMI_result$BMI)),]
}
summary(BMI_result)


Obesity <- read.table("/mnt/nas/gylee/0.GWAS/4.PRS/test/Ukb_obesity_3000_GYQC105/obesity_Case_ids.txt", head=T)
BMI_result$Obesity <- "Control"
wObesity_case <- which(BMI_result$eid %in% Obesity$IID)
BMI_result$Obesity[wObesity_case] <- "Case"

summary(BMI_result[which(BMI_result$Obesity=="Control"),])
summary(BMI_result[which(BMI_result$Obesity=="Case"),])

#
Cli_Obesity_Control <- BMI_result[which(BMI_result$Obesity=="Control"),]
Cli_Obesity_Case <- BMI_result[which(BMI_result$Obesity=="Case"),]

Cli_Obesity_Case <- Cli_Obesity_Case[which(Cli_Obesity_Case$BMI >= 30),]
Cli_Obesity_Case <- Cli_Obesity_Case[which(Cli_Obesity_Case$inch >= 40),]
dim(Cli_Obesity_Case)
summary(Cli_Obesity_Case)
Cli_Obesity_Control <- Cli_Obesity_Control[which(Cli_Obesity_Control$BMI < 24.9),]
Cli_Obesity_Control <- Cli_Obesity_Control[which(Cli_Obesity_Control$inch < 40),]
dim(Cli_Obesity_Control)
summary(Cli_Obesity_Control)

Cli_Obesity <- rbind(Cli_Obesity_Case, Cli_Obesity_Control)
Cli_Obesity$Obesity <- gsub("Control",0,Cli_Obesity$Obesity)
Cli_Obesity$Obesity <- gsub("Case",1,Cli_Obesity$Obesity)
Cli_Obesity$Obesity <- as.numeric(Cli_Obesity$Obesity)

Cli_Obesity_m <- matchit(Obesity ~ Age + Gender, ratio=2, method="nearest", data=Cli_Obesity)
matched_Cli_Obesity <- match.data(Cli_Obesity_m)
matched_Case <- matched_Cli_Obesity[which(matched_Cli_Obesity$Obesity==1),]
matched_Control <- matched_Cli_Obesity[which(matched_Cli_Obesity$Obesity==0),]
summary(matched_Case)
summary(matched_Control)

PSM_Case <- matched_Case[,c("eid","Obesity")]
PSM_Control <- matched_Control[,c("eid","Obesity")]
PSM <- rbind(PSM_Case, PSM_Control)
colnames(PSM) <- c("IID","Pheno")
write.table(PSM,"/mnt/nas/gylee/0.GWAS/Ukb/Obesity_BMI_WC_pheno.txt", row.names=F, col.names=T, quote=F)
writeLines(PSM$IID, "/mnt/nas/gylee/0.GWAS/Ukb/Obesity_BMI_WC.list")
