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
source("/data/keeyoung/gy_RNA/Code/source/source_RNAseq_go.R")
source("/data/keeyoung/gy_RNA/Code/source/source_RNAseq_kegg.R")

print("Make_DIR")
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
gtf <- read.table("/data/keeyoung/gy_RNA/REFERENCE/gtf/edit_Homo_sapiens.GRCh38.100.gtf", sep="\t", head=T)
#go_category
go_cg <- read.table("/data/keeyoung/gy_RNA/REFERENCE/GOterm/go_categories.txt", sep=",", head =F)

#For gene_symbol
o <- order(gtf$gene_name, decreasing=TRUE)
gtf <- gtf[o,]
d <- duplicated(gtf$gene_name)
gtf <- gtf[!d,]
symbol_gtf <- gtf
rownames(symbol_gtf) <- rownames(gtf)
#symbol_gtf <- symbol_gtf[,-1]

symbol_gtf$Length <- as.numeric(symbol_gtf$Length)
symbol_gtf$Start <- as.numeric(symbol_gtf$Start)
symbol_gtf$End <- as.numeric(symbol_gtf$End)
#symbol_function <- symbol_gtf$function.
#symbol_gtf <- symbol_gtf[,-4]

cnt <- iCD
#Matching
symbol_gtf <- na.omit(symbol_gtf[match(rownames(cnt), rownames(symbol_gtf)),])
cnt <- na.omit(cnt[match(rownames(symbol_gtf), rownames(cnt)),])

y <- DGEList(cnt, genes=symbol_gtf)
print({
	paste0("Raw dimension is ", dim(y))
	})
#Ensembl to GeneSymbol
rownames(y) <- y$genes$gene_name

#Ensembl add function and extract proten_coding.
g <- grep("protein_coding", y$genes$function.)
y <- y[g,]

#Unique Gene Symbol
o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
d <- duplicated(y$genes$gene_name)
y <- y[!d,]

#Remove LOC gene
remove <- grep("^LOC", y$genes$gene_name, value=T)
y_uniq <- y[which(!y$genes$gene_name %in% remove),]

print({
	paste0("Extract Protein Codind Genes and Remove LOC genes : ", dim(y_uniq))
	})
#Normalized
keep <- filterByExpr(y_uniq)
y_uniq <- y_uniq[keep, , keep.lib.sizes=FALSE]
y_uniq <- y_uniq[filterByExpr(y_uniq),]
y_uniq <- calcNormFactors(y_uniq)
y_uniq <- estimateDisp(y_uniq,sample.design)
#y_uniq <- estimateDisp(y_uniq) #Using classic mode
print({
	paste0("Do Norm is :" ,dim(y_uniq))
	})

#####Analysis####
#PCA
y_uniq.lcpm <- cpm(y_uniq, log=TRUE)
y_uniq.lcpm.t <- t(y_uniq.lcpm)
y_uniq.pca <- pca(y_uniq.lcpm.t, ncomp=10, center=TRUE, scale=FALSE)
png(filename=paste0(work_dir,"01.Plot/allgenes_pca.png"), width=20, height=15, units="cm", res=200)
pca.variates <- y_uniq.pca$variates[[1]]
pca.variates <- pca.variates[,1:2]
pca.variates <- cbind(pca.variates, sample.pheno)
pca.variates <- as.data.frame(pca.variates)
pca.variates$PC1 <- as.numeric(pca.variates$PC1)
pca.variates$PC2 <- as.numeric(pca.variates$PC2)
p <- ggplot(pca.variates, aes(PC1,PC2, shape=factor(sample.pheno)))
p <- p + geom_point(aes(colour = factor(sample.pheno)), size=5)
p <- p + scale_color_manual(values=c('#56B4E9','#E69F00'))
#p < p + scale_x_continuous(limits = c(-150,150)) + scale_y_continuous(limits = c(-150,150))
p <- p + coord_cartesian(xlim = c(min(pca.variates$PC1)-5, max(pca.variates$PC1)+5))
print({p})
dev.off()

print("Run glmQy_uniqLF")
#glmQLF = F-test
#my.contrasts <- makeContrasts(class=NR-R, levels=sample.design)
Ffit <- glmQLFit(y_uniq, sample.design)
Ftest <- glmQLFTest(Ffit, contrast=my.contrasts[,"class"])
print("Ftest")
print({dim(Ftest)})
png(filename=paste0(work_dir,"01.Plot/glmQLF.png"), width=60, height=50, units="cm", res=200)
print({plotMD(Ftest) + abline(h=c(-1, 1), col="blue")})
dev.off()
write.table(Ftest$table,paste0(work_dir,"02.Table/",name,"_Ftest_table.txt"), quote=F, sep='\t')

#Pathway analysis of Ftest genes
print("Pathway analyzed Ftest")
GO_func(work_dir,Ftest)
KEGG_func(work_dir,Ftest)

idfound <- Ftest$gene$gene_name %in% mappedRkeys(org.Hs.egSYMBOL)
Ftest <- Ftest[idfound,]
egENTREZID <- toTable(org.Hs.egSYMBOL)
entrez_m <- match(Ftest$genes$gene_name , egENTREZID$symbol)
Ftest$genes$entrezID <- egENTREZID$gene_id[entrez_m]

#DecideTests and Plot
dt <- (decideTests(Ftest, p.value=0.05, adjust.method="fdr", lfc=1))
w <- which(dt==0)
qcFtest <- Ftest[-w,]
print({
	paste0("qcFtest's demension is : ", dim(qcFtest))
	})
png(filename=paste0(work_dir,"01.Plot/QC_glmQLF_L.png"), width=60, height=50, units="cm", res=200)
print({plotMD(qcFtest) + abline(h=c(-1, 1), col="blue")})
dev.off()
write.table(qcFtest$table,paste0(work_dir,"02.Table/",name,"_qcFtest_table.txt"), quote=F, sep='\t')

idfound <- qcFtest$gene$gene_name %in% mappedRkeys(org.Hs.egSYMBOL)
qcFtest <- qcFtest[idfound,]
egENTREZID <- toTable(org.Hs.egSYMBOL)
entrez_m <- match(qcFtest$genes$gene_name , egENTREZID$symbol)
qcFtest$genes$entrezID <- egENTREZID$gene_id[entrez_m]
David <- cbind(qcFtest$genes$entrezID ,qcFtest$table$logFC, qcFtest$gene$gene_name)
print({
	paste0("David's dimention is : ", dim(David))
	})
write.table(David,paste0(work_dir,"02.Table/",name,"_qcFtest_table_David.txt"), quote=F, sep=',')

#Genes that pass through fdr and logFC
w <- which(rownames(y_uniq) %in% rownames(qcFtest$genes))
y_uniq_dt <- y_uniq[w,]
y_uniq_dt_heatmap <- y_uniq_dt
rownames(y_uniq_dt_heatmap) <- y_uniq_dt$genes$symbol
#y_uniq_dt.lcpm <- cpm(y_uniq_dt, log=T)
y_uniq_dt.lcpm <- cpm(y_uniq_dt_heatmap, log=T)
print({
	paste0("PCA gene's dimension is ", dim(y_uniq_dt.lcpm))
	})

#PCA
y_uniq_dt.lcpm.t <- t(y_uniq_dt.lcpm)
y_uniq_dt.pca <- pca(y_uniq_dt.lcpm.t, ncomp=10, center=TRUE, scale=FALSE)
png(filename=paste0(work_dir,"01.Plot/dt_pca.png"), width=20, height=15, units="cm", res=200)
pca.variates <- y_uniq_dt.pca$variates[[1]]
pca.variates <- pca.variates[,1:2]
pca.variates <- cbind(pca.variates, sample.pheno)
pca.variates <- as.data.frame(pca.variates)
pca.variates$PC1 <- as.numeric(pca.variates$PC1)
pca.variates$PC2 <- as.numeric(pca.variates$PC2)
p <- ggplot(pca.variates, aes(PC1,PC2, shape=factor(sample.pheno)))
p <- p + geom_point(aes(colour = factor(sample.pheno)), size=5)
p <- p + scale_color_manual(values=c('#56B4E9','#E69F00'))
#p < p + scale_x_continuous(limits = c(-150,150)) + scale_y_continuous(limits = c(-150,150))
p <- p + coord_cartesian(xlim = c(min(pca.variates$PC1)-5, max(pca.variates$PC1)+5))
print({p})
dev.off()
write.table(y_uniq_dt.pca$variates, paste0(work_dir,"01.Plot/pca_variates.csv"), quote=F, col.names=T, row.names=T,sep=",")

#Heatmap
for(i in 1:length(rownames(y_uniq_dt.lcpm))){
y_uniq_dt.lcpm[i,] <- scale(y_uniq_dt.lcpm[i,], center=TRUE, scale=FALSE)
}
cols <- list(pheno=c("CD_in"="blue", "CD_no"="orange"))
anno <- HeatmapAnnotation(
    pheno=sample.pheno, col=cols,
    simple_anno_size = unit(2.5, "cm"),height = unit(2.5, "cm"),
    annotation_name_rot = 45,
    annotation_name_gp=gpar(fontsize=50),
    annotation_legend_param=list(
		pheno=list(title_gp=gpar(fontsize=40), labels_gp=gpar(fontsize=35))
	)
)

png(filename=paste0(work_dir,"01.Plot/Heatmap.png"), width=60, height=150, units="cm", res=200)
print({Heatmap(y_uniq_dt.lcpm, top_annotation=anno,
        column_names_gp=grid::gpar(fontsize=12), row_names_gp = grid::gpar(fontsize = 12),
        width=unit(45, "cm"), height=unit(135,"cm"),
        heatmap_legend_param = list(title = "mat", labels_gp = gpar(fontsize = 25), title_gp = gpar(fontsize = 30))
)})
dev.off()



}

