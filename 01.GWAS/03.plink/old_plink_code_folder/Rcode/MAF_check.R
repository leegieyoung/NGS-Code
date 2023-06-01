maf_freq <- read.table("merge_control_freq.frq", header =TRUE, as.is=T)
pdf("MAF.pdf")
hist(maf_freq[,5],main = "MAF distribution", xlab = "MAF")
dev.off()


