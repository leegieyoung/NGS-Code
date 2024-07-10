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
CD_infla_L <- grep("CD_infla_L",sample.meta$V7)
CD_infla_R <- grep("CD_infla_R",sample.meta$V7)
CD_normal_L <- grep("CD_normal_L",sample.meta$V7)
sample.meta <- sample.meta[c(CD_infla_L, CD_infla_R, CD_normal_L),]
sample.meta$V7 <- gsub("CD_normal_L","CD_uninf_L", sample.meta$V7)
v4 <- grep("Truseq", sample.meta$V4)
sample.meta <- sample.meta[v4,]
sample.ID <- sample.meta$V1
sample.ID <- gsub("-",".", sample.ID)
sample.pheno <- sample.meta$V7 #Only_CD_INvsNO_L_R
#sample.pheno <- sample.meta$V3 #Likely DR,NDR
#sample.pheno <- gsub("^L$","L_R",sample.pheno)
#sample.pheno <- gsub("^R$","L_NR",sample.pheno)

names(sample.pheno) <- sample.ID
sample.group <- factor(sample.pheno)
m <- match(names(sample.pheno),colnames(iCD))
iCD <- iCD[,m]

if(all.equal(names(sample.pheno),colnames(iCD))=="TRUE"){
sample.group <- factor(sample.pheno)
sample.design <- model.matrix(~0 + sample.group)
colnames(sample.design) <- levels(sample.group)
}

#my.contrasts <- makeContrasts(class=L_NR-L_R, levels=sample.design) #L_NR, L_R
my.contrasts <- makeContrasts(class=CD_infla_L-CD_infla_R, levels=sample.design) #CD_in, CD_no

#source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_lcpm.R"))
#source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_lcpm_ensembleToGenesymbol.R"))
#source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_cpm_ensembleToGenesymbol.R"))
#source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_normalized_ensemblToGenesymbol.R")) #220921
#source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_normalized_ensemblToGenesymbol_CD_in_Common.R")) #221003
source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_normalized_ensemblToGenesymbol_analysis_CD_INvsNO_L_R.R")) #221010
RNAseq("CD_INvsNO_Tru_L_R",iCD)

