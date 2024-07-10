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
rownames(symbol_gtf) <- gtf$gene_name
symbol_gtf <- symbol_gtf[,-1]

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

#Ensembl add function and extract proten_coding.
g <- grep("protein_coding", y$genes$function.)
y <- y[g,]

#Unique Gene Symbol
o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
d <- duplicated(y$genes)
y <- y[!d,]

#Remove LOC gene
remove <- grep("^LOC", rownames(y), value=T)
y_uniq <- y[which(!rownames(y) %in% remove),]

#Normalized
#keep <- filterByExpr(y_uniq)
#y_uniq <- y_uniq[keep, , keep.lib.sizes=FALSE]
#y_uniq <- y_uniq[filterByExpr(y_uniq),]
#y_uniq <- calcNormFactors(y_uniq)
#y_uniq <- estimateDisp(y_uniq,sample.design)
#y_uniq <- estimateDisp(y_uniq) #Using classic mode

###
print("End of loading raw data")
print("Do edgeR and Annotation using org.Hs")
y_uniq.cpm <- cpm(y_uniq)
y_uniq.lcpm <- cpm(y_uniq, log=TRUE)
write.table(y_uniq.cpm, paste0(work_dir,"02.Table","/",name,"_cpm.txt"), sep='\t', quote=F,col.names=T, row.names=T)
write.table(y_uniq.lcpm, paste0(work_dir,"02.Table","/",name,"_lcpm.txt"), sep='\t', quote=F,col.names=T, row.names=T)

y_uniq.rpkm <- rpkm(y_uniq)
y_uniq.lrpkm <- rpkm(y_uniq,log=TRUE)
write.table(y_uniq.rpkm, paste0(work_dir,"02.Table","/",name,"_rpkm.txt"), sep='\t', quote=F,col.names=T, row.names=T)
write.table(y_uniq.lrpkm, paste0(work_dir,"02.Table","/",name,"_lrpkm.txt"), sep='\t', quote=F,col.names=T, row.names=T)

}

