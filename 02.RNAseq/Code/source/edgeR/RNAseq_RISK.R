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

#iCD <- read.table("/data/keeyoung/gy_RNA/02.featureCounts/02.hg38_hisat2_result/RISK.txt",sep='\t', head=T)
iCD <- read.table("/data/keeyoung/gy_RNA/02.featureCounts/02.hg38_hisat2_result/RISK_199.txt",sep='\t', head=T)
#colnames(iCD) <- gsub(".txt","",colnames(iCD))

#sample.meta <- read.table("/data/keeyoung/scRNA/iCD/01.Samplelist/RISK_cohort.txt", sep='\t')
sample.meta <- read.table("/data/keeyoung/scRNA/iCD/01.Samplelist/RISK_CD_inflamed.txt", sep='\t')
colnames(sample.meta) <- sample.meta[1,]
sample.meta <- sample.meta[-1,]
rownames(sample.meta) <- sample.meta[,1]
sample.meta <- sample.meta[,-1]

#Extract
#dataset <- grep("RISK",sample.meta$dataset)
sample.meta <- sample.meta[which(sample.meta$dataset %in% "RISK"),]
sample.meta <- sample.meta[which(sample.meta$status %in% "Inflamed"),]
#status <- grep("Inflamed",sample.meta$status)
#NR <- grep("^NR$", sample.meta$response)
#R <- grep("^R$", sample.meta$response)
#sample.meta <- sample.meta[c(R,NR),]

sample.ID <- sample.meta$SRA
sample.pheno <- sample.meta$response
names(sample.pheno) <- sample.ID
sample.pheno[is.na(sample.pheno)] <- "Unkwown"
m <- match(names(sample.pheno),colnames(iCD))
iCD <- iCD[,m]

if(all.equal(names(sample.pheno), colnames(iCD))=="TRUE"){
#sample.group <- factor(sample.pheno, levels=c("NSNP","stricture","penetrating"), order = T) more 3 pheno.
sample.group <- factor(sample.pheno)

#sample.design must have two phenotypes.
sample.design <- model.matrix(~0 + sample.group)
colnames(sample.design) <- levels(sample.group)
}

Cell_frac <- read.table("/data/keeyoung/gy_RNA/01.Samplelist/MGI_RNR_CD_in.txt", sep='\t', head=T)
rownames(Cell_frac) <- Cell_frac[,1]
Cell_frac <- Cell_frac[,-1]
Cell_frac <- Cell_frac[which(rownames(Cell_frac) %in% names(sample.pheno)),]

my.contrasts <- makeContrasts(class=NR-R, levels=sample.design)

#source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_lcpm.R"))
#source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_lcpm_ensembleToGenesymbol.R"))
#source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_cpm_ensembleToGenesymbol.R"))
source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_normalized_ensemblToGenesymbol.R")) #220921
#RNAseq("RISK",iCD)
RNAseq("RISK_199",iCD)

