suppressMessages({
library(edgeR)
library(limma)
library(RColorBrewer)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(repr)
library(statmod)
library(GO.db)
library(ComplexHeatmap)
library(mixOmics)
})

RNAseq <- function(name,iCD){
#DIR
print("Make_DIR")
pre_dir <- "/data/keeyoung/gy_RNA/iCD/"
if (!dir.exists(paste0(pre_dir,name,"/"))){
	dir.create(paste0(pre_dir,name,"/"))
}
work_dir <- paste0(pre_dir,name,"/")

if (!dir.exists(paste0(work_dir,"01.Plot"))){
	dir.create(paste0(work_dir,"01.Plot"))
}
if (!dir.exists(paste0(work_dir,"02.Table"))){
    dir.create(paste0(work_dir,"02.Table"))
}
if (!dir.exists(paste0(work_dir,"03.Pathway"))){
    dir.create(paste0(work_dir,"03.Pathway"))
}

#Loading data
print("Loading data")
#iCD <- read.table("/data/keeyoung/scRNA/iCD/output/convertFeatureCounts/iCD_9patients.txt", sep=' ', head=T)
#gtf <- read.table("/data/keeyoung/gy_RNA/REFERENCE/gtf/anno_Homo_sapiens.GRCh38.100.gtf", sep=" ", head=T)
#go_category
go_cg <- read.table("/data/keeyoung/gy_RNA/REFERENCE/GOterm/go_categories.txt", sep=",", head =F)

rownames(iCD) <- iCD[,1]
iCD <- iCD[,-1]
sample.meta <- read.table("/data/keeyoung/scRNA/iCD/01.Samplelist/sample.csv", sep=',')
sample.ID <- sample.meta$V1
sample.pheno <- sample.meta$V4
sample.pheno <- substr(sample.pheno,1,5)
names(sample.pheno) <- sample.ID
sample.group <- factor(sample.pheno)
sample.group <- factor(sample.pheno)
sample.design <- model.matrix(~0 + sample.group)
colnames(sample.design) <- levels(sample.group)

my.contrasts <- makeContrasts(INvsNO=Infla-Uninf, levels=sample.design)

INvsNO <- c("Infla","Uninf")
wINvsNO <- which(sample.group %in% INvsNO)

col_INvsNO <- list(cINvsNO=c("Infla"="blue","Uninf"="yellow"))

###
print("End of loading raw data")
cnt <- iCD
print("Do edgeR and Annotation using org.Hs")
cnt.lcpm <- cpm(cnt, log=TRUE)
y <- DGEList(cnt, group=sample.group, genes=rownames(cnt.lcpm),)

#Gene Symbol -> Entrezid
idfound <- rownames(y) %in% mappedRkeys(org.Hs.egSYMBOL)
y <- y[idfound,]
egENTREZID <- toTable(org.Hs.egSYMBOL)
m <- match(rownames(y) , egENTREZID$symbol)
y$genes$symbol <- y$genes$genes

y$genes$genes <- egENTREZID$gene_id[m]

#Uniq Entrezid
o <- order(rowSums(y$counts), decreasing=TRUE)
y_uniq <- y[o,]
d <- duplicated(y_uniq$genes$genes)
y_uniq <- y_uniq[!d,]
rownames(y_uniq$counts) <- rownames(y_uniq$genes) <- y_uniq$genes$genes

#Remove LOC gene
remove <- grep("^LOC", y_uniq$genes$symbol, value=T)
y_uniq <- y_uniq[which(!y_uniq$genes$genes %in% remove),]

keep <- filterByExpr(y_uniq)
y_uniq <- y_uniq[keep, , keep.lib.sizes=FALSE]
y_uniq <- y_uniq[filterByExpr(y_uniq),]
y_uniq <- calcNormFactors(y_uniq)
y_uniq <- estimateDisp(y_uniq,sample.design)
y_uniq.lcpm <- cpm(y_uniq, log=TRUE)

y_uniq.zscore <- y_uniq.lcpm
rownames(y_uniq.zscore) <- y_uniq$genes$genes
y_uniq.zscore <- t(scale(t(y_uniq.zscore)))

print("Run glmQLF")
#glmQLF = F-test
Ffit <- glmQLFit(y_uniq, sample.design)
Ftest <- glmQLFTest(Ffit, contrast=my.contrasts[,"INvsNO"])
print("Ftest")
print({dim(Ftest)})
png(filename=paste0(work_dir,"01.Plot/glmQLF.png"), width=60, height=50, units="cm", res=200)
print({plotMD(Ftest) + abline(h=c(-1, 1), col="blue")})
dev.off()

#DecideTests and Plot
dt <- (decideTests(Ftest, p.value=0.05, adjust.method="fdr", lfc=1))
w <- which(dt==0)
qcFtest <- Ftest[-w,]
print("qcFtest")
print({dim(qcFtest)})
png(filename=paste0(work_dir,"01.Plot/dt_glmQLF.png"), width=60, height=50, units="cm", res=200)
print({plotMD(qcFtest) + abline(h=c(-1, 1), col="blue")})
dev.off()

#Genes that pass through fdr and logFC
w <- which(rownames(y_uniq) %in% rownames(qcFtest$genes))
y_uniq_dt <- y_uniq[w,]
y_uniq_dt_heatmap <- y_uniq_dt
rownames(y_uniq_dt_heatmap) <- y_uniq_dt$genes$symbol
#y_uniq_dt.lcpm <- cpm(y_uniq_dt, log=T)
y_uniq_dt.lcpm <- cpm(y_uniq_dt_heatmap, log=T)

print("Run PCA")
#PCA
y_uniq_dt.lcpm.t <- t(y_uniq_dt.lcpm)
y_uniq_dt.pca <- pca(y_uniq_dt.lcpm.t, ncomp=10, center=TRUE, scale=FALSE)
png(filename=paste0(work_dir,"01.Plot/test_pca.png"), width=20, height=15, units="cm", res=200)
pca.variates <- y_uniq_dt.pca$variates[[1]]
pca.variates <- pca.variates[,1:2]
pca.variates <- cbind(pca.variates, sample.pheno)
pca.variates <- as.data.frame(pca.variates)
pca.variates$PC1 <- as.numeric(pca.variates$PC1)
pca.variates$PC2 <- as.numeric(pca.variates$PC2)
p <- ggplot(pca.variates, aes(PC1,PC2, shape=factor(sample.pheno)))
p <- p + geom_point(aes(colour = factor(sample.pheno)), size=5)
#p <- p + scale_color_manual(values=c('#56B4E9','#E69F00'))
p <- p + scale_color_manual(values=c('#56B4E9','#E69F00','#000000','#
EEEEEE'))
#p < p + scale_x_continuous(limits = c(-150,150)) + scale_y_continuous(limits = c(-150,150))
p <- p + coord_cartesian(xlim = c(min(pca.variates$PC1)-5, max(pca.variates$PC1)+5))
print({p})
dev.off()
write.table(y_uniq_dt.pca$variates, paste0(work_dir,"01.Plot/pca_variates.csv"), quote=F, col.names=T, row.names=T,sep=",")


#Heatmap
for(i in 1:length(rownames(y_uniq_dt.lcpm))){
y_uniq_dt.lcpm[i,] <- scale(y_uniq_dt.lcpm[i,], center=TRUE, scale=FALSE)
}
anno <- HeatmapAnnotation(
	INvsNO=sample.pheno, col=col_INvsNO,
	simple_anno_size = unit(2.5, "cm"),height = unit(2.5, "cm"),
	annotation_name_rot = 45,
	annotation_name_gp=gpar(fontsize=50),
	annotation_legend_param=list(INvsNO=list(title_gp=gpar(fontsize=40), labels_gp=gpar(fontsize=35)))
	)
png(filename=paste0(work_dir,"01.Plot/Heatmap.png"), width=60, height=150, units="cm", res=200)
Heatmap(y_uniq_dt.lcpm, top_annotation=anno,
        column_names_gp=grid::gpar(fontsize=12), row_names_gp = grid::gpar(fontsize = 12),
        width=unit(45, "cm"), height=unit(135,"cm"),
        heatmap_legend_param = list(title = "mat", labels_gp = gpar(fontsize = 25), title_gp = gpar(fontsize = 30))
)
dev.off() 


#Online-tools
table <- cbind(qcFtest$table, qcFtest$genes)
write.table(table, paste0(work_dir,"02.Table/genes.csv"), quote=F, col.names=T, row.names=T,sep=",")
GOnet <- cbind(qcFtest$genes$genes ,qcFtest$table$logFC)
write.table(GOnet, paste0(work_dir,"02.Table/GOnet.csv"), quote=F, col.names=T, row.names=T,sep=",")

David <- cbind(rownames(qcFtest$genes) ,qcFtest$table$logFC)
write.table(David, paste0(work_dir,"02.Table/David.csv"), quote=F, col.names=T, row.names=T,sep=",")


#Pathway-kegg
keg.all <- kegga(qcFtest, species="Hs")
cutoff.keg <- subset(keg.all, keg.all$P.Up < 0.05 | keg.all$P.Down < 0.05)
cutoff.keg$Pathway <- gsub(",", "|", cutoff.keg$Pathway)
cutoff.keg$Pathway <- gsub(" ", "_", cutoff.keg$Pathway)
keg.up.down <- cbind(rownames(cutoff.keg), cutoff.keg$Pathway, cutoff.keg$P.Up, cutoff.keg$P.Down)
write.table(keg.up.down, paste0(work_dir,"03.Pathway/","kegg.up.down.txt"), quote=F, col.names=F, row.names=F, sep="\t")

#Pathway-GO
go.all <- goana(qcFtest, species="Hs")
cutoff.go <- subset(go.all, go.all$P.Up < 0.05 | go.all$P.Down < 0.05)
cutoff.go$Term <- gsub(" ", "_", cutoff.go$Term)
cutoff.go$Term <- gsub(",", "|", cutoff.go$Term)
go.up.down <- cbind(rownames(cutoff.go), cutoff.go$Term, cutoff.go$P.Up, cutoff.go$P.Down)
write.table(keg.up.down, paste0(work_dir,"03.Pathway/","go.up.down.txt"), quote=F, col.names=F, row.names=F, sep="\t")

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
write.table(go.up.bp, paste0(work_dir,"03.Pathway/","go.up.bp.txt"),quote=F, col.names=F, row.names=F, sep="\t")
go.up.cc <- go.up[which(go.up[,1] %in% cc[,1]),]
write.table(go.up.cc,  paste0(work_dir,"03.Pathway/","go.up.cc.txt"),quote=F, col.names=F, row.names=F, sep="\t")
go.up.mf <- go.up[which(go.up[,1] %in% mf[,1]),]
write.table(go.up.mf, paste0(work_dir,"03.Pathway/","go.up.mf.txt"),quote=F, col.names=F, row.names=F, sep="\t")

go.down.bp <- go.down[which(go.down[,1] %in% bp[,1]),]
write.table(go.down.bp, paste0(work_dir,"03.Pathway/","go.down.bp.txt"),quote=F, col.names=F, row.names=F, sep="\t")
go.down.cc <- go.down[which(go.down[,1] %in% cc[,1]),]
write.table(go.down.cc, paste0(work_dir,"03.Pathway/","go.down.cc.txt"),quote=F, col.names=F, row.names=F, sep="\t")
go.down.mf <- go.down[which(go.down[,1] %in% mf[,1]),]
write.table(go.down.mf, paste0(work_dir,"03.Pathway/","go.down.mf.txt"),quote=F, col.names=F, row.names=F, sep="\t")

}

