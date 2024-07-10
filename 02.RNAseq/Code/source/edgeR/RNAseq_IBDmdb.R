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

IBDmdb <- read.delim("/data/keeyoung/gy_RNA/IBDmdb/host_tx_counts.tsv", header=T)
IBDmdb <- IBDmdb[,-which(colnames(IBDmdb) %in% c("CSMDRVXI","MSM719ME"))] #Zero size
pre_dir <- "/data/keeyoung/gy_RNA/06.output/"
name="IBDmdb"

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
if (!dir.exists(paste0(work_dir,"04.WGCNA"))){
    dir.create(paste0(work_dir,"04.WGCNA"))
}

#Loading GTF
print("Loading data")
gtf <- read.table("/data/keeyoung/gy_RNA/REFERENCE/gtf/edit_Homo_sapiens.GRCh38.100.gtf", sep="\t", head=T)

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

cnt <- IBDmdb
#Matching
symbol_gtf <- na.omit(symbol_gtf[match(rownames(cnt), symbol_gtf$gene_name),])
cnt <- na.omit(cnt[match(symbol_gtf$gene_name, rownames(cnt)),])

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


#Normalized
keep <- filterByExpr(y_uniq)
y_uniq <- y_uniq[keep, , keep.lib.sizes=FALSE]
y_uniq <- y_uniq[filterByExpr(y_uniq),]
y_uniq <- calcNormFactors(y_uniq)
y_uniq <- estimateDisp(y_uniq) #Using classic mode

print({
    paste0("Do Norm is :" ,dim(y_uniq))
    })
y_uniq_cpm <- cpm(y_uniq, log=F)
write.table(y_uniq_cpm, paste0(work_dir,"02.Table/",name,"_cpm.txt"), quote=F, sep='\t')

