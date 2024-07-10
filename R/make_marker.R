SAMPLE <- Sys.getenv("SAMPLE")
RESULT_DIR <- Sys.getenv("RESULT_DIR")
CUTOFF <- Sys.getenv("CUTOFF")

prune <- read.table(paste0(RESULT_DIR,"prune/",CUTOFF,".glm"), head=T, comment.char = "")
prune$OR <- as.numeric(prune$OR)
prune$beta <- log(prune$OR)
prune$absBETA <- abs(prune$beta)
#ABS 보정
prune$absEA <- ""
for (A in 1:dim(prune)[1]){
        if(prune$beta[A]>=0){
                prune$absEA[A] <- prune$A1[A]
        } else if(prune$beta[A] < 0){
                if(prune$A1[A]==prune$ALT[A]){
                        prune$absEA[A] <- prune$REF[A]
                } else if(prune$A1[A]!=prune$ALT[A]){
                        prune$absEA[A] <- prune$ALT[A]}
}
}
marker <- prune[,c("ID","X.CHROM","POS","absBETA","absEA")]
write.table(marker, paste0(RESULT_DIR,SAMPLE,".",CUTOFF,".prune_marker.txt"), quote=F, row.names=F, sep='\t', col.names=F)


