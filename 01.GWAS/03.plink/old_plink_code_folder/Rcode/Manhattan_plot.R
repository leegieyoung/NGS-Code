#install.packages("qqman",repos="https://ftp.harukasan.org/CRAN/",lib="/home/keeyoung/Rpackage" ) # location of installation can be changed but has to correspond with the library location 
library("qqman",lib.loc="/home/keeyoung/Rpackage")  
#results_log <- read.table("logistic_results.assoc_2.logistic", head=TRUE)
#jpeg("Logistic_manhattan.jpeg")
#manhattan(results_log,chr="CHR",bp="BP",p="P",snp="SNP", main = "Manhattan plot: logistic")
#dev.off()

#results_as <- read.table("MAF_merge_CCSNP_merge_control_case_NoNA_assoc.assoc", head=TRUE)
#pdf("assoc_manhattan.pdf")
#manhattan(results_as,chr="CHR",bp="BP",p="P",snp="SNP", main = "Manhattan plot", ylim = c(0, 50), cex = 0.6, cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F, chrlabs = c(1:22))
#dev.off()  

#results_as <- read.table("low_pval_merge_CCSNP_merge_control_case_NoNA_assoc.assoc", head=TRUE)
#pdf("assoc_manhattan_annotation.pdf")
#manhattan(results_as,chr="CHR",bp="BP",p="P",snp="SNP", main = "Manhattan plot_annotation", ylim = c(0, 60), cex = 0.6, cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F, chrlabs = c(1:22), annotatePval = 0.00000000619, annotateTop = FALSE)
#dev.off()


#e55_noheader_remove_wrong_MaMi_MAF_merge_CCSNP_merge_control_case_NoNA_assoc.assoc
#results_as <- read.table("e55_noheader_remove_wrong_MaMi_MAF_merge_CCSNP_merge_control_case_NoNA_assoc.assoc", head=TRUE)
#pdf("e55_assoc_manhattan_annotation.pdf")
#manhattan(results_as,chr="CHR",bp="BP",p="P",snp="SNP", main = "Manhattan plot_annotation", ylim = c(0, 60), cex = 0.6, cex.axis = 0.9, col = c("blue4", "orange3"), suggestiveline = F, genomewideline = F, chrlabs = c(1:22), annotatePval = 0.00000000619, annotateTop = FALSE)
#dev.off()

results_as <- read.table("Change", head=TRUE)
png(filename="Manhattan_assoc.png")
manhattan(results_as,chr="CHR",bp="BP",p="P",snp="SNP", ylim = c(0, 10), cex = 0.6, cex.axis = 0.9, col = c("blue4", "orange3"), genomewideline = -log10(5.63e-9) , chrlabs = c(1:22), suggestiveline = F, annotatePval = 0.00001, annotateTop = FALSE)
dev.off()

results_as <- read.table("Change", head=TRUE)
png(filename="QQ-Plot_assoc.png")
qq(results_as$P, col = "blue4", ylim=c(0, 10))
dev.off()
