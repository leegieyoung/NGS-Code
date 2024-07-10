KEGG_func <- function(work_dir,Ftest){
print("Start KEGG_func analysis")
print(work_dir)

#Pathway-GO
idfound <- Ftest$gene$gene_name %in% mappedRkeys(org.Hs.egSYMBOL)
Ftest <- Ftest[idfound,]
egENTREZID <- toTable(org.Hs.egSYMBOL)
m <- match(Ftest$gene$gene_name , egENTREZID$symbol)

Ftest$genes$entrezID <- egENTREZID$gene_id[m]
rownames(Ftest) <- Ftest$genes$entrezID

keg.all <- kegga(Ftest, species="Hs")
cutoff.keg <- subset(keg.all, keg.all$P.Up < 0.05 | keg.all$P.Down < 0.05)
cutoff.keg$Pathway <- gsub(",", "|", cutoff.keg$Pathway)
cutoff.keg$Pathway <- gsub(" ", "_", cutoff.keg$Pathway)
keg.up.down <- cbind(rownames(cutoff.keg), cutoff.keg$Pathway, cutoff.keg$P.Up, cutoff.keg$P.Down)
colnames(keg.up.down) <- c("KEGG","Pathway","Up","Down")
write.table(keg.up.down, paste0(work_dir,"03.Pathway/","kegg.up.down.txt"), quote=F, sep="\t")

}
