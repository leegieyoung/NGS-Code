suppressMessages({
library(edgeR)
library(limma)
library(RColorBrewer)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)
library(repr)
library(statmod)
library(GO.db)
library(ComplexHeatmap)
library(mixOmics)
library(dplyr)
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

y <- DGEList(cnt, genes=symbol_gtf, group=sample.platform)
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
y_uniq_cpm_NoQC <- cpm(y_uniq, log=F)
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
y_uniq_cpm <- cpm(y_uniq, log=F)
write.table(y_uniq_cpm, paste0(work_dir,"02.Table/",name,"_cpm.txt"), quote=F, sep='\t')

#PCA
y_uniq.lcpm <- cpm(y_uniq, log=TRUE)
y_uniq.lcpm.t <- t(y_uniq.lcpm)
y_uniq.pca <- pca(y_uniq.lcpm.t, ncomp=10, center=TRUE, scale=FALSE)
png(filename=paste0(work_dir,"01.Plot/allgenes_pca_No_group.png"), width=20, height=15, units="cm", res=200)
pca.variates <- y_uniq.pca$variates[[1]]
pca.variates <- pca.variates[,1:2]
pca.variates <- cbind(pca.variates, sample.pheno)
pca.variates <- as.data.frame(pca.variates)
pca.variates$PC1 <- as.numeric(pca.variates$PC1)
pca.variates$PC2 <- as.numeric(pca.variates$PC2)
p <- ggplot(pca.variates, aes(PC1,PC2, shape=factor(sample.pheno)))
p <- p + geom_point(aes(colour = factor(sample.pheno)), size=5)
p <- p + scale_color_manual(values=c('#56B4E9','#E69F00',"#AB5757","#82C787"))
#p < p + scale_x_continuous(limits = c(-150,150)) + scale_y_continuous(limits = c(-150,150))
p <- p + coord_cartesian(xlim = c(min(pca.variates$PC1)-5, max(pca.variates$PC1)+5))
print({p})
dev.off()

}

