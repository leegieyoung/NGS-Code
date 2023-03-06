#--
library(edgeR)
library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)

iCD <- read.table("/data/keeyoung/gy_RNA/02.featureCounts/CRC/CRC.txt", sep='\t', head=T)
gtf <- read.table("/data/keeyoung/gy_RNA/REFERENCE/gtf/anno_Homo_sapiens.GRCh38.100.gtf", sep=" ", head=T)
rownames(iCD) <- iCD[,1]
iCD <- iCD[,-1]
cnt <- iCD
cnt.lcpm <- cpm(cnt, log=TRUE)
y <- cnt.lcpm
y <- as.data.frame(y)
idfound <- rownames(y) %in% mappedRkeys(org.Hs.egENSEMBL)
y <- y[idfound,]
dim(y)
m <- match(rownames(y), gtf$ensembl_id)
egENTREZID <- toTable(org.Hs.egENSEMBL)
m <- match(rownames(y), egENTREZID$ensembl_id)
idfound <- egENTREZID$gene_id[m]
y$genes <- ''
y$genes <- idfound
y_order <- y[order(y$genes, decreasing=TRUE),]

d <- duplicated(y_order$genes)
y_uniq <- y_order[!d,]
dim(y_uniq)
egSYMBOL <- toTable(org.Hs.egSYMBOL)
m <- match(y_uniq$genes, egSYMBOL$gene_id)
y_uniq$symbol <- ''
y_uniq$symbol <- egSYMBOL$symbol[m]
remove <- grep("^LOC",y_uniq$symbol, value=T)
y_uniq1 <- y_uniq[which(!y_uniq$symbol %in% remove),]
y_uniq2 <- y_uniq1
rownames(y_uniq2) <- y_uniq1$symbol
y_uniq2 <- y_uniq2[,1:73]
write.table(y_uniq2, "/data/keeyoung/gy_RNA/CRC/CRC_lcpm_GeneSymbol.txt", sep='\t', quote=F, row.names=T, col.names=T)

