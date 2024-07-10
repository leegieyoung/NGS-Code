SAMPLE <- Sys.getenv("SAMPLE")
RESULT_DIR <- Sys.getenv("RESULT_DIR")
print({SAMPLE})
print({RESULT_DIR})
VEPoutput <- read.table(paste0(RESULT_DIR,SAMPLE,".snplist.vep_web.txt"),header=T, comment.char="")
FDRresult <- read.table(paste0(RESULT_DIR,SAMPLE,".PHENO1.glm.logistic.hybrid.FDR"),header=T, comment.char="")
PRUNEresult <- read.table(paste0(RESULT_DIR,"prune/","prune.glm"),header=T, comment.char="")

m <- match(VEPoutput$X.Uploaded_variation, FDRresult$ID)	
VEPoutput$P <- FDRresult$P[m]
VEPoutput$OR <- FDRresult$OR[m]
VEPoutput$CHR <- FDRresult$X.CHROM[m]
VEPoutput$POS <- FDRresult$POS[m]
VEPoutput$ALT <- FDRresult$ALT[m] 
VEPoutput$REF <- FDRresult$REF[m]
VEPoutput$A1 <- FDRresult$A1[m]

VEPoutput$Gene <- gsub("-","None",VEPoutput$Gene)
gene.list <- unique(sort(VEPoutput$Gene))
rmw <- which(gene.list %in% "None")
gene.list <- gene.list[-rmw]

#Total
FDRresult$beta <- log(FDRresult$OR)
total <- FDRresult
total <- total[,c("ID","X.CHROM","POS","A1","ALT","REF", "beta","P")]
colnames(total) <- c("rsID","CHR","POS","A1","ALT","REF","beta","P")

#ABS 보정
total$absEA <- ""
for (A in 1:dim(total)[1]){
        if(total$beta[A]>=0){
                total$absEA[A] <- total$A1[A]
        } else if(total$beta[A] < 0){
                if(total$A1[A]==total$ALT[A]){
                        total$absEA[A] <- total$REF[A]
                } else if(total$A1[A]!=total$ALT[A]){
                        total$absEA[A] <- total$ALT[A]}
}
}
total$absbeta <- abs(total$beta)
write.table(total, paste0(RESULT_DIR,SAMPLE,".total.rmdupID_UniqGene.txt"), quote=F, row.names=F, sep='\t')
marker <- total[,c("rsID","CHR","POS","absbeta","absEA")]
write.table(marker, paste0(RESULT_DIR,SAMPLE,".total_marker.txt"), quote=F, row.names=F, sep='\t', col.names=F)

#Gene
result <- data.frame(col1 = character(), col2 = character(), col3 = numeric(), stringsAsFactors = FALSE)

for (A in 1:length(gene.list)){
	w <- which(VEPoutput$Gene %in% gene.list[A])
	temp <- VEPoutput[w,]
	minP <- temp[which.min(temp$P),]
	minP <- cbind(minP, gene.list[A])
	result <- rbind(result, minP[, c("X.Uploaded_variation","CHR","POS","A1","ALT","REF","gene.list[A]", "OR","P")])
}

d <- which(duplicated(result$X.Uploaded_variation))
if(length(d)>0){
        result <- result[-d,]}

result$OR <- log(result$OR)
colnames(result) <- c("rsID","CHR","POS","A1","ALT","REF","Gene","beta","P")

#ABS 보정
result$absEA <- ""
for (A in 1:dim(result)[1]){
	if(result$beta[A]>=0){
		result$absEA[A] <- result$A1[A]
	} else if(result$beta[A] < 0){
		if(result$A1[A]==result$ALT[A]){
			result$absEA[A] <- result$REF[A]
		} else if(result$A1[A]!=result$ALT[A]){
			result$absEA[A] <- result$ALT[A]}
}
}
result$absbeta <- abs(result$beta)
write.table(result, paste0(RESULT_DIR,SAMPLE,".rmdupID_UniqGene.txt"), quote=F, row.names=F, sep='\t') 

marker <- result[,c("rsID","CHR","POS","absbeta","absEA")]
write.table(marker, paste0(RESULT_DIR,SAMPLE,"_marker.txt"), quote=F, row.names=F, sep='\t', col.names=F)
#Nonew <- which(VEPoutput$Gene %in% "None")
#temp2 <- VEPoutput[Nonew,]
#mapIDw <- which(temp2$X.Uploaded_variation %in% result$rsID)
#temp2 <- temp2[-mapIDw,]

#PRUNE
print({"PRUNE"})
Prumem <- match(PRUNEresult$ID, VEPoutput$X.Uploaded_variation)
VEPoutput <- VEPoutput[Prumem,]
print({"m"})
print({head(m)})

VEPoutput$Gene <- gsub("-","None",VEPoutput$Gene)
gene.list <- unique(sort(VEPoutput$Gene))
rmw <- which(gene.list %in% "None")
gene.list <- gene.list[-rmw]
print({gene.list})

result <- data.frame(col1 = character(), col2 = character(), col3 = numeric(), stringsAsFactors = FALSE)

for (A in 1:length(gene.list)){
        w <- which(VEPoutput$Gene %in% gene.list[A])
        temp <- VEPoutput[w,]
        minP <- temp[which.min(temp$P),]
        minP <- cbind(minP, gene.list[A])
        result <- rbind(result, minP[, c("X.Uploaded_variation","CHR","POS","A1","ALT","REF","gene.list[A]", "OR","P")])
}

d <- which(duplicated(result$X.Uploaded_variation))
if(length(d)>0){
	result <- result[-d,]}

result$OR <- log(result$OR)
colnames(result) <- c("rsID","CHR","POS","A1","ALT","REF","Gene","beta","P")


#ABS 보정
result$absEA <- ""
for (A in 1:dim(result)[1]){
        if(result$beta[A]>=0){
                result$absEA[A] <- result$A1[A]
        } else if(result$beta[A] < 0){
                if(result$A1[A]==result$ALT[A]){
                        result$absEA[A] <- result$REF[A]
                } else if(result$A1[A]!=result$ALT[A]){
                        result$absEA[A] <- result$ALT[A]}
}
}
result$absbeta <- abs(result$beta)
write.table(result, paste0(RESULT_DIR,SAMPLE,".rmdupID_UniqGene.prune.txt"), quote=F, row.names=F, sep='\t')

marker <- result[,c("rsID","CHR","POS","absbeta","absEA")]
write.table(marker, paste0(RESULT_DIR,SAMPLE,".prune_marker.txt"), quote=F, row.names=F, sep='\t', col.names=F)

