GO_func <- function(work_dir,Ftest){
print("Start GO_func analysis")
print(work_dir)
#go_category
go_cg <- read.table("/data/keeyoung/gy_RNA/REFERENCE/GOterm/go_categories.txt", sep=",", head =F)

#Pathway-GO
idfound <- Ftest$gene$gene_name %in% mappedRkeys(org.Hs.egSYMBOL)
Ftest <- Ftest[idfound,]
egENTREZID <- toTable(org.Hs.egSYMBOL)
m <- match(Ftest$gene$gene_name , egENTREZID$symbol)

Ftest$genes$entrezID <- egENTREZID$gene_id[m]
rownames(Ftest) <- Ftest$genes$entrezID

go.all <- goana(Ftest, species="Hs")
cutoff.go <- subset(go.all, go.all$P.Up < 0.05 | go.all$P.Down < 0.05)
cutoff.go$Term <- gsub(" ", "_", cutoff.go$Term)
cutoff.go$Term <- gsub(",", "|", cutoff.go$Term)
#go.up.down <- cbind(rownames(cutoff.go), cutoff.go$Term, cutoff.go$P.Up, cutoff.go$P.Down)
#colnames(go.up.down) <- c("GO","term","Up","Down")
write.table(cutoff.go, paste0(work_dir,"03.Pathway/","go.up.down.txt"), quote=F, sep="\t")

##GO_category
go_namespace <- c("go","namespace")
colnames(go_cg) <- go_namespace
w <- which(go_cg$namespace %in% "biological_process")
bp <- go_cg[w,]
w <- which(go_cg$namespace %in% "molecular_function")
mf <- go_cg[w,]
w <- which(go_cg$namespace %in% "cellular_component")
cc <- go_cg[w,]
go_p <- c("GO","P")

go.up <- cutoff.go[which(cutoff.go$P.Up < 0.05),]
go.up <- cbind(rownames(go.up), go.up$Term, go.up$P.Up)

go.down <- cutoff.go[which(cutoff.go$P.Down < 0.05),]
go.down <- cbind(rownames(go.down), go.down$Term, go.down$P.Down)

go.up.bp <- go.up[which(go.up[,1] %in% bp[,1]),]
write.table(go.up.bp, paste0(work_dir,"03.Pathway/","go.up.bp.txt"),quote=F, sep="\t")
go.up.cc <- go.up[which(go.up[,1] %in% cc[,1]),]
write.table(go.up.cc,  paste0(work_dir,"03.Pathway/","go.up.cc.txt"),quote=F, sep="\t")
go.up.mf <- go.up[which(go.up[,1] %in% mf[,1]),]
write.table(go.up.mf, paste0(work_dir,"03.Pathway/","go.up.mf.txt"),quote=F, sep="\t")

go.down.bp <- go.down[which(go.down[,1] %in% bp[,1]),]
write.table(go.down.bp, paste0(work_dir,"03.Pathway/","go.down.bp.txt"),quote=F, sep="\t")
go.down.cc <- go.down[which(go.down[,1] %in% cc[,1]),]
write.table(go.down.cc, paste0(work_dir,"03.Pathway/","go.down.cc.txt"),quote=F, sep="\t")
go.down.mf <- go.down[which(go.down[,1] %in% mf[,1]),]
write.table(go.down.mf, paste0(work_dir,"03.Pathway/","go.down.mf.txt"),quote=F, sep="\t")
}
