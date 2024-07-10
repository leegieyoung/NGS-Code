#install.packages("qqman",repos="http://cran.cnr.berkeley.edu/",lib="~" ) # location of installation can be changed but has to correspond with the library location 
library("qqman") 
input <- "NR_PRCR_result.assoc.logistic_manhattan_convert"
results_log <- read.table(paste0("/src/data/",input,".txt"), head=TRUE)
png(paste0("/src/data/",input,".png"), res=300, width=3000, height=1800)
manhattan(results_log,chr="CHR",bp="BP",p="P",snp="SNP", main = "Manhattan plot: logistic", 
	  col=c("blue4","orange3"), suggestiveline=-log10(1e-5), genomewideline= -log10(5e-8),
	  chrlabs=c(1:22), annotatePval=1e-5, annotateTop=FALSE)
dev.off()




