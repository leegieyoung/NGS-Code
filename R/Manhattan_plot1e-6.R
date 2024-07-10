#install.packages("qqman",repos="http://cran.cnr.berkeley.edu/",lib="~" ) # location of installation can be changed but has to correspond with the library location 
library("qqman") 
SAMPLE <- Sys.getenv("SAMPLE")
RESULT_DIR <- Sys.getenv("RESULT_DIR")
#var <- scan(paste0(RESULT_DIR,"variant.info"), what=character(0))
CUTOFF <- as.numeric(1e-6)
L_CUTOFF <- -log10(CUTOFF)

results_log <- read.table(paste0(RESULT_DIR,SAMPLE,".PHENO1.glm.logistic.hybrid"), head=TRUE, comment.char = "")
#test <- results_log
print({dim(results_log)})
#rm NA
if(length(which(is.na(results_log$P))) >0){
	results_log <- results_log[-which(is.na(results_log$P)),]}

results_log$P <- as.numeric(results_log$P)
png(paste0(RESULT_DIR,SAMPLE,".png"), res=300, width=3000, height=1800)
manhattan(results_log,chr="X.CHROM",bp="POS",p="P",snp="ID", main = "Manhattan plot: logistic", ylim = c(0, abs(min(log(results_log$P,10)))+1 ), 
	  col=c("blue4","orange3"), suggestiveline=-log10(1e-5), genomewideline= L_CUTOFF,
	  chrlabs=c(1:22), annotatePval=CUTOFF, annotateTop=TRUE)
dev.off()

w <- which(results_log[,"P"] < CUTOFF)
result_FDR <- results_log[w,]
write.table(result_FDR, paste0(RESULT_DIR,SAMPLE,".PHENO1.glm.logistic.hybrid.FDR"), col.names=T, row.names=F, quote=F)



