library(ggplot2);
SAMPLE <- Sys.getenv("SAMPLE")
QC_DIR <- Sys.getenv("QC_DIR")

print({QC_DIR})
print({SAMPLE})

maf_freq <- read.table(paste0(QC_DIR,SAMPLE,"_g_m_maf_hwe_freq.afreq"), header=TRUE, as.is=T)
pdf(paste0(QC_DIR,"MAF_distribution.pdf"))
p <- ggplot(maf_freq, aes(x=ALT_FREQS)) + geom_histogram(binwidth=0.01)
#hist(maf_freq[,6],main = "MAF distribution", xlab = "MAF", breaks=0.01)
dev.off()


