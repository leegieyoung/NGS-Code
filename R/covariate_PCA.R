#!/bin/Rscript
library(ggplot2);
SAMPLE <- Sys.getenv("SAMPLE")
RESULT_DIR <- Sys.getenv("RESULT_DIR")

print({"---"})
print({"PCA"})
print({"---"})

PC20 <- read.table(paste0(RESULT_DIR,SAMPLE,"_PCA.eigenvec"), header=T, comment.char = "")
head(PC20)
print({dim(PC20)})
COVAR <- PC20
colnames(COVAR) <- gsub("X.IID","#IID",colnames(PC20))
write.table(COVAR, paste0(RESULT_DIR,"covariate_forplink_PCA.txt"), row.names=F, col.names=T, sep='\t', quote=F)

