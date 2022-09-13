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
v2 <- grep("CD_infla",sample.meta$V2)
sample.meta <- sample.meta[v2,]
v4 <- grep("Totalseq", sample.meta$V4)
sample.meta <- sample.meta[v4,]

sample.ID <- sample.meta$V1
sample.ID <- gsub("-",".", sample.ID)
sample.pheno <- sample.meta$V2
sample.pheno <- substr(sample.pheno,1,5)
names(sample.pheno) <- sample.ID

#Unknown <- grep("Unknown",sample.pheno)
#sample.pheno <- sample.pheno[-Unknown]

#sample.group <- factor(sample.pheno, levels=c("NSNP","stricture","penetrating"), order = T)
sample.group <- factor(sample.pheno)

#sample.design must have two phenotypes.
#sample.design <- model.matrix(~0 + sample.group)
#colnames(sample.design) <- levels(sample.group)

m <- match(names(sample.pheno),colnames(iCD))
iCD <- iCD[,m]

#source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_lcpm.R"))
source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_lcpm_ensembleToGenesymbol.R"))
RNAseq("CD_in",iCD)
