#install.packages("qqman",repos="http://cran.cnr.berkeley.edu/",lib="~" ) # location of installation can be changed but has to correspond with the library location 
library("qqman") 
library("data.table")
SAMPLE <- Sys.getenv("SAMPLE")
RESULT_DIR <- Sys.getenv("RESULT_DIR")
CUTOFF <- as.numeric("5e-8")
L_CUTOFF <- -log10(CUTOFF)

result <- fread(paste0(RESULT_DIR,"merge.result"), head = TRUE, sep=' ')

colnames(result)[which(colnames(result)=="p.value")] <- "P"
colnames(result)[which(colnames(result)=="CHR")] <- "X.CHROM"
colnames(result)[which(colnames(result)=="MarkerID")] <- "ID"

results_log <- result
#test <- results_log
print({dim(results_log)})
#rm NA
if(length(which(is.na(results_log$P))) >0){
	results_log <- results_log[-which(is.na(results_log$P)),]}
results_log$Total <- (results_log$N_case + results_log$N_ctrl)

write.table(results_log, paste0(RESULT_DIR,SAMPLE,".result"), col.names=T, row.names=F, quote=F)


w0 <- which(results_log$P==0)
if(length(w0)!=0){
	results_log$P[w0] <- 0.05
}

results_log$P <- as.numeric(results_log$P)

print({head(results_log)})
png(paste0(RESULT_DIR,SAMPLE,".png"), res=300, width=3000, height=1800)
manhattan(results_log,chr="X.CHROM",bp="POS",p="P",snp="ID", main = "Manhattan plot: SAIGE", ylim = c(0, abs(round(min(log(results_log$P,10)),10))+1 ), 
	  col=c("blue4","orange3"), suggestiveline=-log10(1e-5), genomewideline= L_CUTOFF,
	  chrlabs=c(1:22), annotatePval=CUTOFF, annotateTop=TRUE)
dev.off()

#
if(length(w0)!=0){
        results_log$P[w0] <- 1e-999
}

w <- which(results_log[,"P"] < CUTOFF)
result_FDR <- results_log[w,]
write.table(result_FDR, paste0(RESULT_DIR,SAMPLE,".PHENO1.glm.logistic.hybrid.FDR"), col.names=T, row.names=F, quote=F)

#1e-5
CUTOFF <- as.numeric(1e-5)
L_CUTOFF <- -log10(CUTOFF)

print({CUTOFF})
print({L_CUTOFF})

w <- which(results_log[,"P"] < CUTOFF)
result_suggest <- results_log[w,]
write.table(result_suggest, paste0(RESULT_DIR,SAMPLE,".PHENO1.glm.logistic.hybrid.suggest"), col.names=T, row.names=F, quote=F)


#1e-4
CUTOFF <- as.numeric(1e-4)
L_CUTOFF <- -log10(CUTOFF)

print({CUTOFF})
print({L_CUTOFF})

w <- which(results_log[,"P"] < CUTOFF)
result_suggest <- results_log[w,]
write.table(result_suggest, paste0(RESULT_DIR,SAMPLE,".PHENO1.glm.logistic.hybrid.1e-4"), col.names=T, row.names=F, quote=F)

