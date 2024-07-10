SAMPLE <- Sys.getenv("SAMPLE")
QC_DIR <- Sys.getenv("QC_DIR")

het <- read.table(paste0(QC_DIR,SAMPLE,"_R_check.het"), header=T, comment.char = "")
pdf(paste0(QC_DIR,"heterozygosity.pdf"))
het$HET_RATE = (het$OBS_CT - het$O.HOM.)/het$OBS_CT
hist(het$HET_RATE, xlab="Heterozygosity Rate", ylab="Frequency", main= "Heterozygosity Rate")
dev.off()
