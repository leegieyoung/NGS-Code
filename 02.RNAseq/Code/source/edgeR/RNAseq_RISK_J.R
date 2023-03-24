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
name="RISK_J"
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

#
#RISK <- read.table("/data/keeyoung/gy_RNA/GIMAT/GSE134881/GSE134881_rmdupli_LOC_RISKfpkmAll_forR.txt",sep='\t')
Genes <- data.frame(t(RISK[which(rownames(RISK) %in% c("TNF","IL6","FAP")),]))
Genes <- log(Genes,2)
colnames(Genes) <- c("logFAP","logTNF","logIL6")
#
sample.meta <- read.table("/data/keeyoung/scRNA/iCD/01.Samplelist/199_RISK_CD_inflamed.csv", sep=',')
colnames(sample.meta) <- sample.meta[1,]
sample.meta <- sample.meta[-1,]
rownames(sample.meta) <- sample.meta[,1]
sample.meta <- sample.meta[,-1]
sample.meta$response[is.na(sample.meta$response)] <- "Unknwon"
Genes$response <- sample.meta$response[which(rownames(sample.meta) %in% rownames(Genes))]

#Plot
#rm Overexpression, J10300
Genes <- Genes[-which(rownames(Genes) %in% "J10300"),]

png(filename=paste0(work_dir,"01.Plot/TNF_IL6.png"), width=60, height=50, units="cm", res=200)
ggplot(data=Genes, mapping=aes(x=logTNF,y=logIL6, shape=response, color=response)) + 
	geom_point(size=5) + 
	geom_smooth(method=lm, se=FALSE) + 
	xlab(expression("log"[2]*"TNF")) +
	ylab(expression("log"[2]*"IL6")) +
	theme(axis.title=element_text(size=25), axis.text=element_text(size=20), legend.title=element_text(size=25), legend.text=element_text(size=20)) +
	guides(colour= guide_legend(override.aex=list(size=1.5)))
dev.off()

png(filename=paste0(work_dir,"01.Plot/FAP_IL6.png"), width=60, height=50, units="cm", res=200)
ggplot(Genes, aes(x=logFAP,y=logIL6, shape=response, color=response)) + 
    geom_point(size=5) +
    geom_smooth(method=lm, se=FALSE) +
    xlab(expression("log"[2]*"FAP")) +
    ylab(expression("log"[2]*"IL6")) +
    theme(axis.title=element_text(size=25), axis.text=element_text(size=20), legend.title=element_text(size=25), legend.text=element_text(size=20)) +
    guides(colour= guide_legend(override.aex=list(size=1.5)))
dev.off()
