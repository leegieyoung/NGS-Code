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
gtf <- read.table("/data/keeyoung/gy_RNA/REFERENCE/gtf/anno_Homo_sapiens.GRCh38.100.gtf", sep=" ", head=T)
#go_category
go_cg <- read.table("/data/keeyoung/gy_RNA/REFERENCE/GOterm/go_categories.txt", sep=",", head =F)

###
print("End of loading raw data")
cnt <- iCD
print("Do edgeR and Annotation using org.Hs")
cnt.lcpm <- cpm(cnt, log=TRUE)
y <- DGEList(cnt, group=sample.group, genes=rownames(cnt.lcpm),)

idfound <- rownames(y) %in% mappedRkeys(org.Hs.egENSEMBL)
y <- y[idfound,]

#Ensembl add function and extract proten_coding.
m <- match(rownames(y), gtf$ensembl_id)
Func <- gtf$function.[m]
y$genes <- cbind(y$genes, Func)
g <- grep("protein_coding", y$genes$Func)
y <- y[g,]

egENTREZID <- toTable(org.Hs.egENSEMBL)
m <- match(rownames(y), egENTREZID$ensembl_id)
y$genes$entrezID <- egENTREZID$gene_id[m]

#Entrezid -> Gene Symbol
egSYMBOL <- toTable(org.Hs.egSYMBOL)
m <- match(y$genes$entrezID, egSYMBOL$gene_id)
y$genes$symbol <- egSYMBOL$symbol[m]

#Remove LOC gene
remove <- grep("^LOC", y$genes$symbol, value=T)
y <- y[which(!y$genes$symbol %in% remove),]

#Remove duplicated GeneSymbol
o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
d <- duplicated(y$genes$symbol)
y_uniq <- y[!d,]
rownames(y_uniq) <- y_uniq$genes$symbol

print("Do edgeR and Annotation using org.Hs")
cnt.cpm <- cpm(y_uniq, log=TRUE)
write.table(cnt.cpm, paste0(work_dir,"02.Table","/",name,"_lcpm_Genesymbol_QC.txt"), sep='\t', quote=F,col.names=T, row.names=T)

}

