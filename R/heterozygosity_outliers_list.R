SAMPLE <- Sys.getenv("SAMPLE")
QC_DIR <- Sys.getenv("QC_DIR")

het <- read.table(paste0(QC_DIR,SAMPLE,"_R_check.het"), header=T, comment.char = "")
het$HET_RATE = (het$OBS_CT - het$O.HOM.)/het$OBS_CT
het_fail = subset(het, (het$HET_RATE < mean(het$HET_RATE)-3*sd(het$HET_RATE)) | (het$HET_RATE > mean(het$HET_RATE)+3*sd(het$HET_RATE)));
het_fail$HET_DST = (het_fail$HET_RATE-mean(het$HET_RATE))/sd(het$HET_RATE);
write.table(het_fail, paste0(QC_DIR,SAMPLE,"_fail-het-qc.txt"), row.names=FALSE)
