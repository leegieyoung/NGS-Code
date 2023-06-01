hwe<-read.table (file="merge_control_hwe.hwe", header=TRUE)
pdf("histhwe.pdf")
hist(hwe[,9],main="Histogram HWE")
dev.off()

