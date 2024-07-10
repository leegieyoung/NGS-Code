install.packages("qqman",repos="https://ftp.harukasan.org/CRAN/",lib="~" ) # location of installation can be changed but has to correspond with the library location 
library("qqman",lib.loc="~") 
#results_log <- read.table("logistic_results.7_assoc_2.logistic", head=TRUE)
#jpeg("QQ-Plot_logistic.jpeg")
#qq(results_log$P, main = "Q-Q plot of GWAS p-values : log")
#dev.off()

results_as <- read.table("7_assoc_results.assoc", head=TRUE)
pdf("QQ-Plot_7_assoc.pdf")
qq(results_as$P, main = "Q-Q plot of GWAS p-values : log")
dev.off()

