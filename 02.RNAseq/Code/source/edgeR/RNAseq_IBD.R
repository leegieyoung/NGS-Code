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

pre_dir <- "/data/keeyoung/gy_RNA/06.output/"

iCD <- read.table("/data/keeyoung/gy_RNA/02.featureCounts/gene_id/454.IBD.fc_P_M.merge",sep='\t', head=T)
colnames(iCD) <- gsub(".TR","",colnames(iCD))

sample.meta <- read.table("/data/keeyoung/scRNA/iCD/01.Samplelist/sample.csv", sep=',')
#Extract
sample.ID <- sample.meta$V1
sample.ID <- gsub("-",".", sample.ID)
#sample.pheno <- sample.meta$V2 #Only_CD_in
sample.pheno <- sample.meta$V2 #CD_in, CD_no, UC_in, UC_no

sample.pheno <- substr(sample.pheno,1,5)
names(sample.pheno) <- sample.ID

    #platform
sample.platform <- sample.meta$V4
names(sample.platform) <- sample.meta$V1
names(sample.platform) <- gsub("-",".", names(sample.platform))
sample.platform <- sample.platform[which(names(sample.platform) %in% names(sample.pheno))]
sample.platform <- factor(sample.platform)

sample.group <- factor(sample.pheno)
m <- match(names(sample.pheno),colnames(iCD))
iCD <- iCD[,m]

if(all.equal(names(sample.pheno),colnames(iCD))=="TRUE"){
print("All equal")
sample.group <- factor(sample.pheno)
sample.design <- model.matrix(~0 + sample.group)
colnames(sample.design) <- levels(sample.group)
}

#my.contrasts <- makeContrasts(class=L_NR-L_R, levels=sample.design)

#source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_normalized_ensemblToGenesymbol_CD_in_Article.R")) ##221011
source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_normalized_ensemblToGenesymbol_IBD.R")) ##221011
#source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_normalized_ensemblToGenesymbol_IBD_filter_Group.R")) #221108
RNAseq("IBD",iCD)


