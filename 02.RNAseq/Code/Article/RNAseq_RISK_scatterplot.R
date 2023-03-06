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

#my.contrasts <- makeContrasts(class=NR-R, levels=sample.design)
print("Make_DIR")
name="RISK_199"
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
gtf <- read.table("/data/keeyoung/gy_RNA/REFERENCE/gtf/edit_Homo_sapiens.GRCh38.100.gtf", sep="\t", head=T)
#go_category
go_cg <- read.table("/data/keeyoung/gy_RNA/REFERENCE/GOterm/go_categories.txt", sep=",", head =F)

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
#symbol_function <- symbol_gtf$function.
#symbol_gtf <- symbol_gtf[,-4]

cnt <- iCD
#Matching
symbol_gtf <- na.omit(symbol_gtf[match(rownames(cnt), rownames(symbol_gtf)),])
cnt <- na.omit(cnt[match(rownames(symbol_gtf), rownames(cnt)),])

y <- DGEList(cnt, genes=symbol_gtf)

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
#y_uniq <- estimateDisp(y_uniq,sample.design)
y_uniq <- estimateDisp(y_uniq) #Using classic mode

#220930_for RISK_cohort, scatter plot
RISK <- cpm(y_uniq, log=TRUE)

Genes <- data.frame(t(RISK[which(rownames(RISK) %in% c("THBS2","TNF","FAP")),]))
colnames(Genes) <- c("logTHBS2","logTNF","logFAP")
Genes$response <- sample.pheno[which(names(sample.pheno) %in% rownames(Genes))]

#

#The scatter plot
png(filename=paste0(work_dir,"01.Plot/TNF_FAP.png"), width=60, height=50, units="cm", res=200)
ggplot(data=Genes, mapping=aes(x=logTNF,y=logFAP, shape=response, color=response)) +
    geom_point(size=5) +
    geom_smooth(method=lm, se=FALSE) +
    xlab(expression("log"[2]*"TNF")) +
    ylab(expression("log"[2]*"FAP")) +
    scale_x_continuous(limits = c(min(Genes$logTNF)-1, max(Genes$logTNF)+1))+
    scale_y_continuous(limits = c(min(Genes$logFAP)-1, max(Genes$logFAP)+1)) +
    theme(axis.title=element_text(size=25), axis.text=element_text(size=20), legend.title=element_text(size=25), legend.text=element_text(size=20)) +
    guides(colour= guide_legend(override.aex=list(size=1.5)))
dev.off()

png(filename=paste0(work_dir,"01.Plot/THBS2_FAP.png"), width=60, height=50, units="cm", res=200)
ggplot(Genes, aes(x=logTHBS2,y=logFAP, shape=response, color=response)) +
    geom_point(size=5) +
    geom_smooth(method=lm, se=FALSE) +
    xlab(expression("log"[2]*"THBS2")) +
    ylab(expression("log"[2]*"FAP")) +
    scale_x_continuous(limits = c(min(Genes$logTHBS2)-1, max(Genes$logTHBS2)+1))+
    scale_y_continuous(limits = c(min(Genes$logFAP)-1, max(Genes$logFAP)+1)) +
    theme(axis.title=element_text(size=25), axis.text=element_text(size=20), legend.title=element_text(size=25), legend.text=element_text(size=20)) +
    guides(colour= guide_legend(override.aex=list(size=1.5)))
dev.off()

RNR_Genes <- Genes[which(Genes$response %in% c("R","NR")),]
#scatter plot
png(filename=paste0(work_dir,"01.Plot/RNR_TNF_FAP.png"), width=60, height=50, units="cm", res=200)
ggplot(data=RNR_Genes, mapping=aes(x=logTNF,y=logFAP, shape=response, color=response)) +
    geom_point(size=5) +
    geom_smooth(method=lm, se=FALSE) +
    xlab(expression("log"[2]*"TNF")) +
    ylab(expression("log"[2]*"FAP")) +
    scale_x_continuous(limits = c(min(Genes$logTNF)-1, max(Genes$logTNF)+1))+
    scale_y_continuous(limits = c(min(Genes$logFAP)-1, max(Genes$logFAP)+1)) +
    theme(axis.title=element_text(size=25), axis.text=element_text(size=20), legend.title=element_text(size=25), legend.text=element_text(size=20)) +
    guides(colour= guide_legend(override.aex=list(size=1.5)))
dev.off()

png(filename=paste0(work_dir,"01.Plot/RNR_THBS2_FAP.png"), width=60, height=50, units="cm", res=200)
ggplot(RNR_Genes, aes(x=logTHBS2,y=logFAP, shape=response, color=response)) +
    geom_point(size=5) +
    geom_smooth(method=lm, se=FALSE) +
    xlab(expression("log"[2]*"THBS2")) +
    ylab(expression("log"[2]*"FAP")) +
    scale_x_continuous(limits = c(min(Genes$logTHBS2)-1, max(Genes$logTHBS2)+1))+
    scale_y_continuous(limits = c(min(Genes$logFAP)-1, max(Genes$logFAP)+1)) +
    theme(axis.title=element_text(size=25), axis.text=element_text(size=20), legend.title=element_text(size=25), legend.text=element_text(size=20)) +
    guides(colour= guide_legend(override.aex=list(size=1.5)))
dev.off()

