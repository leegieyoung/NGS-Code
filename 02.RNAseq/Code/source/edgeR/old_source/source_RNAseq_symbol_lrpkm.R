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

RNAseq_rpkm <- function(name,iCD){
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
#iCD <- read.table("/data/keeyoung/scRNA/iCD/output/sctype/cibersortx/iCD_count.txt", sep='\t', head=T)
gtf <- read.table("/data/keeyoung/gy_RNA/REFERENCE/gtf/edit_Homo_sapiens.GRCh38.100.gtf", sep="\t", head=T)
#go_category
go_cg <- read.table("/data/keeyoung/gy_RNA/REFERENCE/GOterm/go_categories.txt", sep=",", head =F)

#For gene_symbol
o <- order(gtf[,1], decreasing=TRUE)
gtf <- gtf[o,]
d <- duplicated(gtf[,1])
gtf <- gtf[!d,]
symbol_gtf <- as.data.frame(cbind(gtf[,1], gtf[,4]))
colnames(symbol_gtf) <- c("gene_symbol","length")
symbol_gtf[,2] <- as.numeric(symbol_gtf[,2])

cnt <- iCD
#Matching
symbol_gtf <- na.omit(symbol_gtf[match(rownames(cnt), symbol_gtf[,1]),])
cnt <- na.omit(cnt[match(symbol_gtf[,1], rownames(cnt)),])

y <- DGEList(cnt, genes=symbol_gtf)

#Ensembl add function and extract proten_coding.
m <- match(rownames(y), gtf$gene_name)
Func <- gtf$function.[m]
y$genes <- cbind(y$genes, Func)
g <- grep("protein_coding", y$genes$Func)
y <- y[g,]

#Unique Gene Symbol
o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
d <- duplicated(y$genes$gene_symbol)
y <- y[!d,]

#Remove LOC gene
remove <- grep("^LOC", rownames(y), value=T)
y_uniq <- y[which(!rownames(y) %in% remove),]

y_uniq.rpkm <- rpkm(y_uniq, log=TRUE)
write.table(y_uniq.rpkm, paste0(work_dir,"02.Table/","iCD_lfpkm_matchedGTF.txt"), quote=F, col.names=T, row.names=T, sep="\t")

#
#Remove worth-wihole genes when processing filterByExpr
#keep <- filterByExpr(y_uniq)
#y_uniq <- y_uniq[keep, , keep.lib.sizes=FALSE]
#y_uniq <- y_uniq[filterByExpr(y_uniq),]
#y_uniq <- calcNormFactors(y_uniq)

#y_uniq <- estimateDisp(y_uniq)

#y_uniq_rm_worth-gene.rpkm <- rpkm(y_uniq)

#write.table(y_uniq_rm_worth-gene.rpkm, paste0(work_dir,"02.Table/","rm_worth-gene_fpkm.txt"), quote=F, col.names=T, row.names=F, sep="\t")

}
