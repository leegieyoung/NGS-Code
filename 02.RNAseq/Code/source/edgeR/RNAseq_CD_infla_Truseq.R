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
v4 <- grep("Truseq", sample.meta$V4)
sample.meta <- sample.meta[v4,]

sample.ID <- sample.meta$V1
sample.ID <- gsub("-",".", sample.ID)
#sample.pheno <- sample.meta$V2 #Only_CD_in
sample.pheno <- sample.meta$V3 #Likely DR,NDR
sample.pheno <- gsub("^L$","L_R",sample.pheno)
sample.pheno <- gsub("^R$","L_NR",sample.pheno)

sample.pheno <- substr(sample.pheno,1,5)
names(sample.pheno) <- sample.ID

#Unknown <- grep("Unknown",sample.pheno)
#sample.pheno <- sample.pheno[-Unknown]

#sample.group <- factor(sample.pheno, levels=c("NSNP","stricture","penetrating"), order = T)
sample.group <- factor(sample.pheno)

m <- match(names(sample.pheno),colnames(iCD))
iCD <- iCD[,m]

if(all.equal(names(sample.pheno),colnames(iCD))=="TRUE"){
sample.group <- factor(sample.pheno)
sample.design <- model.matrix(~0 + sample.group)
colnames(sample.design) <- levels(sample.group)
}

my.contrasts <- makeContrasts(class=L_NR-L_R, levels=sample.design)

#source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_lcpm.R"))
#source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_lcpm_ensembleToGenesymbol.R"))
#source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_cpm_ensembleToGenesymbol.R"))
#source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_normalized_ensemblToGenesymbol.R")) #220921
#source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_normalized_ensemblToGenesymbol_CD_in_Common.R")) #221003
#RNAseq("CD_in_Tru",iCD)
#source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_normalized_ensemblToGenesymbol_CD_in_Tru_rmCD_INvsNO.R")) ##221006 GO func
#RNAseq("CD_in_Tru",iCD)
#source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_normalized_ensemblToGenesymbol_CD_in_Tru_rmCD_INvsNO_After_221011.R")) #2210111 After norm, Common genes
#RNAseq("CD_in_Tru_After",iCD)

#source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_normalized_ensemblToGenesymbol_CD_in_Article.R"))
#RNAseq("CD_in_Tru_article",iCD)

#source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_normalized_ensemblToGenesymbol_CD_in_rmCD_INvsNO_702genes.R"))
#RNAseq("CD_in_Tru_792genes",iCD)
source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_normalized_ensemblToGenesymbol_CD_in_Article_246genes.R")) ##221130
RNAseq("CD_in_Tru_article",iCD)
