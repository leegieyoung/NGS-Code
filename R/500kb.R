library("stringr")

SAMPLE <- Sys.getenv("SAMPLE")
RESULT_DIR <- Sys.getenv("RESULT_DIR")
CUTOFF <- Sys.getenv("CUTOFF")
#SAIGE_DIR <- Sys.getenv("SAIGE_DIR")
SAIGE_DIR <- "/mnt/nas/gylee/0.GWAS/2.SAIGE_result/"
result <- read.table(paste0(SAIGE_DIR, "Ukb_Train_Obesity/result/Ukb_Train_Obesity.PHENO1.glm.logistic.hybrid.FDR"), head=T)

colnames(result)[which(colnames(result) %in% "X.CHROM")] <- "CHR"
colnames(result)[which(colnames(result) %in% "Allele1")] <- "A1"
colnames(result)[which(colnames(result) %in% "Allele2")] <- "A2"
colnames(result)[which(colnames(result) %in% "P")] <- "P.value"

#Only SNP and Indel
wREMOVE <- ""
for (A in 1:length(rownames(result))){
        if(abs(nchar(result$A1[A])-nchar(result$A2[A])) > 1){
        print({abs(nchar(result$A1[A])-nchar(result$A2[A]))})
        wREMOVE <- c(wREMOVE,A)
}
}
wREMOVE <- as.numeric(wREMOVE[-1])
result <- result[-wREMOVE,]

win_result <- as.data.frame(matrix(ncol=length(colnames(result)),nrow=0))
colnames(win_result) <- colnames(result)

#minP in window 500kb
for (i in 1:22){
w <- which(result$CHR==i)
 
chr_result <- result[w,]
 
for (A in 1:100){
if(length(rownames(chr_result))==0){
        break
} else {
#Add SNP of min P
win_result <- rbind(win_result,chr_result[which.min(chr_result$P.value),])
 
#Remove window
MIN_POS <- as.numeric(chr_result[which.min(chr_result$P.value),]$POS)
if(MIN_POS - 50000 < 0){
        L_MIN_POS <-  0
} else{
        L_MIN_POS <- MIN_POS - 50000
}
R_MIN_POS <- MIN_POS + 50000
 
chr_result <- subset(chr_result, POS <= L_MIN_POS | POS >= R_MIN_POS)
}
}
}
dim(win_result)

write.table(win_result, paste0(SAIGE_DIR, "Ukb_Train_Obesity/result/Ukb_Train_Obesity.PHENO1.glm.logistic.hybrid.FDR.win50.minP"), row.names=F, col.names=T, quote=F)
writeLines(win_result$ID, paste0(SAIGE_DIR, "Ukb_Train_Obesity/result/Ukb_Train_Obesity.PHENO1.glm.logistic.hybrid.FDR.win50.minP.IID"))

