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

QC_sPLSDA <- function(name,iCD){
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


#my.contrasts <- makeContrasts(INvsNO=Infla-Uninf, levels=sample.design)

#INvsNO <- c("Infla","Uninf")
#wINvsNO <- which(sample.group %in% INvsNO)

#col_INvsNO <- list(cINvsNO=c("Infla"="blue","Uninf"="yellow"))

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

#Remove Duplicated Genesymbols
o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
y_order <- y[order(y$genes$symbol, decreasing=TRUE),]
d <- duplicated(y_order$genes$symbol)
y_uniq <- y[!d,]

#edgeR filter
keep <- filterByExpr(y_uniq)
y_uniq <- y_uniq[keep, , keep.lib.sizes=FALSE]
y_uniq <- y_uniq[filterByExpr(y_uniq),]
y_uniq <- calcNormFactors(y_uniq)
y_uniq <- estimateDisp(y_uniq,sample.design)
print({"End of edgeR filtering"})

#CPM
y_uniq.lcpm <- cpm(y_uniq, log=TRUE)
#write.table(y_uniq.lcpm, paste0(work_dir,"02.Table/cpm.csv"), sep='\t', col.names=T, row.names=T, quote=F)

#sPLSDA
y_uniq.lcpm.t <- t(y_uniq.lcpm)

optimal.ncomp <- 1 
optimal.keepX <- 100
y_uniq.final.splsda <- splsda(y_uniq.lcpm.t, sample.group,
						ncomp = optimal.ncomp,
						keepX = optimal.keepX)

print({"Make plot of tune"})
sample.col <- sample.group
sample.col <- gsub("other","#F68B33",sample.col) #orange
sample.col <- gsub("NSNP","#388ECC",sample.col) #blue

#legend=list(legend=levels(sample.group),
legend=list(legend=c("other","NSNP"),
            col=c("#F68B33","#388ECC"),
            title="Last behavior",
            cex=1.4)

png(filename=paste0(work_dir,"01.Plot/sPLSDA_final_heatmap_comp1.png"), width=70, height=30, units="cm", res=200)
print({cim(y_uniq.final.splsda, row.sideColors = sample.col,
			keysize=c(0.2,0.9),
			legend = legend)
	})
dev.off()

#Loading genes
zero <- (which(y_uniq.final.splsda$loadings$X=="0"))
col1 <- rownames(y_uniq.final.splsda$loadings$X)[-zero]
col2 <- as.numeric(y_uniq.final.splsda$loadings$X[-zero])
y_uniq.gene <- as.data.frame(cbind(col1,col2))
y_uniq.gene$col2 <- as.numeric(y_uniq.gene$col2)
y_uniq.gene <- y_uniq.gene[order(y_uniq.gene[,2], decreasing=TRUE),]
write.table(y_uniq.gene, paste0(work_dir,"02.Table/y_uniq.gene.txt"),sep='\t', quote=F)


}

