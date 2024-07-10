#install.packages("qqman",repos="http://cran.cnr.berkeley.edu/",lib="~" ) # location of installation can be changed but has to correspond with the library location 
library("qqman") 
library("data.table")
SAMPLE <- Sys.getenv("SAMPLE")
RESULT_DIR <- Sys.getenv("RESULT_DIR")
var <- scan(paste0(RESULT_DIR,"variant.info"), what=character(0))
#CUTOFF <- as.numeric(0.05/length(var))
CUTOFF <- 5e-8
CUTOFF <- as.numeric(CUTOFF)
L_CUTOFF <- -log10(CUTOFF)

#results_log <- read.table(paste0(RESULT_DIR,SAMPLE,".PHENO1.glm.logistic.hybrid"), head=TRUE, comment.char = "")
#results_log <- read.table(paste0(RESULT_DIR,SAMPLE,".PHENO1.glm.logistic.hybrid"), head=TRUE, comment.char = "",  colClasses = c("numeric","numeric","character","character","character","character","character","character","numeric","character","character","numeric","numeric","numeric","numeric","numeric","character"))
results_log <- fread(paste0(RESULT_DIR,SAMPLE,".PHENO1.glm.logistic.hybrid"), head=TRUE, sep='\t')
#test <- read.table(paste0(RESULT_DIR,"test.txt"), head=TRUE, comment.char = "",  colClasses = c("character","character","character","character","character","character","character","character","numeric","character","character","numeric","numeric","numeric","numeric","numeric","character"))
results_log2 <- results_log 
#test <- results_log
print({dim(results_log)})
#rm NA
if(length(which(is.na(results_log$P))) >0){
	results_log <- results_log[-which(is.na(results_log$P)),]}

w0 <- which(results_log$P==0)
if(length(w0)!=0){
	results_log$P[w0] <- 0.05
}


print({head(results_log)})
png(paste0(RESULT_DIR,SAMPLE,".png"), res=300, width=3000, height=1800)
manhattan(results_log,chr="#CHROM",bp="POS",p="P",snp="ID", main = "Manhattan plot: logistic", ylim = c(0, abs(round(min(log(results_log$P,10)),10))+1 ), 
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

