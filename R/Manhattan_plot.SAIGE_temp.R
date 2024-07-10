#install.packages("qqman",repos="http://cran.cnr.berkeley.edu/",lib="~" ) # location of installation can be changed but has to correspond with the library location 
library("qqman") 
SAMPLE <- Sys.getenv("SAMPLE")
RESULT_DIR <- Sys.getenv("RESULT_DIR")
CUTOFF <- as.numeric("5e-8")
L_CUTOFF <- -log10(CUTOFF)

result <- read.table(paste0("/mnt/nas/gylee/0.GWAS/2.SAIGE_result/", SAMPLE, "/result/merge.result"), head = TRUE)

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
