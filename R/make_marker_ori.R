SAMPLE <- Sys.getenv("SAMPLE")
RESULT_DIR <- Sys.getenv("RESULT_DIR")
prune <- read.table(paste0(RESULT_DIR,"prune/","prune.glm"), head=T, comment.char = "")
prune$OR <- as.numeric(prune$OR)
prune$beta <- log(prune$OR)
prune$absbeta <- abs(prune$beta)
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
marker <- prune[,c("ID","X.CHROM","POS","absbeta","absEA")]
write.table(marker, paste0(RESULT_DIR,SAMPLE,".prune_marker.txt"), quote=F, row.names=F, sep='\t', col.names=F)

min0.02.prune <- read.table(paste0(RESULT_DIR,"prune/","min0.02.prune.glm"), head=T, comment.char = "")
min0.02.prune$OR <- as.numeric(min0.02.prune$OR)
min0.02.prune$beta <- log(min0.02.prune$OR)
min0.02.prune$absbeta <- abs(min0.02.prune$beta)
#ABS 보정
min0.02.prune$absEA <- ""
for (A in 1:dim(min0.02.prune)[1]){
        if(min0.02.prune$beta[A]>=0){
                min0.02.prune$absEA[A] <- min0.02.prune$A1[A]
        } else if(min0.02.prune$beta[A] < 0){
                if(min0.02.prune$A1[A]==min0.02.prune$ALT[A]){
                        min0.02.prune$absEA[A] <- min0.02.prune$REF[A]
                } else if(min0.02.prune$A1[A]!=min0.02.prune$ALT[A]){
                        min0.02.prune$absEA[A] <- min0.02.prune$ALT[A]}
}
}
min0.02.marker <- min0.02.prune[,c("ID","X.CHROM","POS","absbeta","absEA")]
write.table(min0.02.marker, paste0(RESULT_DIR,SAMPLE,".min0.02.prune_marker.txt"), quote=F, row.names=F, sep='\t', col.names=F)


#1e-5
suggest.prune <- read.table(paste0(RESULT_DIR,"prune/","prune.suggest.glm"), head=T, comment.char = "")
suggest.prune$OR <- as.numeric(suggest.prune$OR)
suggest.prune$beta <- log(suggest.prune$OR)
suggest.prune$absbeta <- abs(suggest.prune$beta)
#ABS 보정
suggest.prune$absEA <- ""
for (A in 1:dim(suggest.prune)[1]){
        if(suggest.prune$beta[A]>=0){
                suggest.prune$absEA[A] <- suggest.prune$A1[A]
        } else if(suggest.prune$beta[A] < 0){
                if(suggest.prune$A1[A]==suggest.prune$ALT[A]){
                        suggest.prune$absEA[A] <- suggest.prune$REF[A]
                } else if(suggest.prune$A1[A]!=suggest.prune$ALT[A]){
                        suggest.prune$absEA[A] <- suggest.prune$ALT[A]}
}
}
suggest.marker <- suggest.prune[,c("ID","X.CHROM","POS","absbeta","absEA")]
write.table(suggest.marker, paste0(RESULT_DIR,SAMPLE,".suggest.prune_marker.txt"), quote=F, row.names=F, sep='\t', col.names=F)

#1e-4
suggest2.prune <- read.table(paste0(RESULT_DIR,"prune/","prune.1e-4.glm"), head=T, comment.char = "")
suggest2.prune$OR <- as.numeric(suggest2.prune$OR)
suggest2.prune$beta <- log(suggest2.prune$OR)
suggest2.prune$absbeta <- abs(suggest2.prune$beta)
#ABS 보정
suggest2.prune$absEA <- ""
for (A in 1:dim(suggest2.prune)[1]){
        if(suggest2.prune$beta[A]>=0){
                suggest2.prune$absEA[A] <- suggest2.prune$A1[A]
        } else if(suggest2.prune$beta[A] < 0){
                if(suggest2.prune$A1[A]==suggest2.prune$ALT[A]){
                        suggest2.prune$absEA[A] <- suggest2.prune$REF[A]
                } else if(suggest2.prune$A1[A]!=suggest2.prune$ALT[A]){
                        suggest2.prune$absEA[A] <- suggest2.prune$ALT[A]}
}
}
suggest2.marker <- suggest2.prune[,c("ID","X.CHROM","POS","absbeta","absEA")]
write.table(suggest2.marker, paste0(RESULT_DIR,SAMPLE,".1e-4.prune_marker.txt"), quote=F, row.names=F, sep='\t', col.names=F)

