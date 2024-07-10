SAMPLE <- Sys.getenv("SAMPLE")
RESULT_DIR <- Sys.getenv("RESULT_DIR")
CUTOFF <- Sys.getenv("CUTOFF")
SAIGE_DIR <- Sys.getenv("SAIGE_DIR")

prune <- read.table(paste0(RESULT_DIR,"prune/",CUTOFF,".glm"), head=T, comment.char = "")
prune$absBETA <- abs(prune$BETA)
#ABS 보정
prune$absEA <- ""
for (A in 1:dim(prune)[1]){
        if(prune$BETA[A]>=0){
                prune$absEA[A] <- prune$Allele2[A]
        } else if(prune$BETA[A] < 0){
		prune$absEA[A] <- prune$Allele1[A]
}
}
marker <- prune[,c("ID","X.CHROM","POS","absBETA","absEA")]
head(marker)
write.table(marker, paste0(RESULT_DIR,SAMPLE,".",CUTOFF,".prune_marker.txt"), quote=F, row.names=F, sep='\t', col.names=F)
