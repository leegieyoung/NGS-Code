install.packages("qqman",repos="https://ftp.harukasan.org/CRAN/",lib="~" ) # location of installation can be changed but has to correspond with the library location 
library("qqman",lib.loc="~")  
#results_log <- read.table("logistic_results.7_assoc_2.logistic", head=TRUE)
#jpeg("Logistic_manhattan.jpeg")
#manhattan(results_log,chr="CHR",bp="BP",p="P",snp="SNP", main = "Manhattan plot: logistic")
#dev.off()

results_as <- read.table("7_assoc_results.assoc", head=TRUE)
pdf("7_assoc_manhattan.pdf")
manhattan(results_as,chr="CHR",bp="BP",p="P",snp="SNP", main = "Manhattan plot: 7_assoc")
dev.off()  




