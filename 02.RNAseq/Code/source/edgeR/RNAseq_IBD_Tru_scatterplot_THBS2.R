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
name="IBD_Tru"
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
iCD <- read.table("/data/keeyoung/gy_RNA/02.featureCounts/gene_id/454.IBD.fc_P_M.merge",sep='\t', head=T)
colnames(iCD) <- gsub(".TR","",colnames(iCD))

sample.meta <- read.table("/data/keeyoung/scRNA/iCD/01.Samplelist/sample.csv", sep=',')
#Extract
sample.meta <- sample.meta[which(sample.meta$V4 %in% "Truseq"),]
sample.ID <- sample.meta$V1
sample.ID <- gsub("-",".", sample.ID)
#sample.pheno <- sample.meta$V2 #Only_CD_in
sample.pheno <- sample.meta$V2 #CD_in, CD_no, UC_in, UC_no
sample.pheno <- substr(sample.pheno,1,5)
names(sample.pheno) <- sample.ID

#sample.RNR <- sample.meta$V5 #DR,NDR,Not_CD_in,Unknown
#names(sample.RNR) <- sample.ID
#sample.RNR <- gsub("Not_CD_in","Unknown",sample.RNR)
sample.RNR <- sample.meta$V7 #"CD_infla_L"  "CD_infla_R"  "CD_normal_L" "UC_infla_L"  "UC_normal_L"
names(sample.RNR) <- sample.ID

sample.platform <- sample.meta$V4
names(sample.platform) <- sample.meta$V1
names(sample.platform) <- gsub("-",".", names(sample.platform))
sample.platform <- sample.platform[which(names(sample.platform) %in% names(sample.pheno))]

sample.group <- factor(sample.pheno)
m <- match(names(sample.pheno),colnames(iCD))
iCD <- iCD[,m]

if(all.equal(names(sample.pheno),colnames(iCD))=="TRUE"){
print("All equal")
sample.group <- factor(sample.pheno)
sample.design <- model.matrix(~0 + sample.group)
colnames(sample.design) <- levels(sample.group)
}

#
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

y <- DGEList(cnt, genes=symbol_gtf, group=sample.platform)
print({
    paste0("Raw dimension is ", dim(y))
    })
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
y_uniq_cpm_NoQC <- cpm(y_uniq, log=F)
print({
    paste0("Extract Protein Codind Genes and Remove LOC genes : ", dim(y_uniq))
    })
#Normalized
#keep <- filterByExpr(y_uniq)
#y_uniq <- y_uniq[keep, , keep.lib.sizes=FALSE]
#y_uniq <- y_uniq[filterByExpr(y_uniq),]
y_uniq <- calcNormFactors(y_uniq)
y_uniq <- estimateDisp(y_uniq,sample.design)

y_uniq.lcpm <- cpm(y_uniq, log=TRUE)
Genes <- data.frame(t(y_uniq.lcpm[which(rownames(y_uniq.lcpm) %in% c("THBS2","IL6","FAP")),]))
colnames(Genes) <- c("logFAP","logTHBS2","logIL6")
#
Genes$response <- sample.pheno[which(names(sample.pheno) %in% rownames(Genes))]
Genes$RNR <- sample.RNR[which(names(sample.RNR) %in% rownames(Genes))]
OnlyRNR_Genes <- Genes[which(Genes$RNR %in% c("CD_infla_L","CD_infla_R")),]
OnlyRNR_Genes$RNR <- gsub("CD_infla_L","DR",OnlyRNR_Genes$RNR)
OnlyRNR_Genes$RNR <- gsub("CD_infla_R","NDR",OnlyRNR_Genes$RNR)
#Plot

png(filename=paste0(work_dir,"01.Plot/THBS2_FAP.png"), width=60, height=50, units="cm", res=200)
ggplot(data=Genes, mapping=aes(x=logTHBS2,y=logFAP, shape=response, color=response)) + 
	geom_point(size=5) + 
	geom_smooth(method=lm, se=FALSE) + 
	xlab(expression("log"[2]*"THBS2")) +
	ylab(expression("log"[2]*"FAP")) +
    scale_x_continuous(limits = c(min(Genes$logFAP)-1, max(Genes$logFAP)+1))+
    scale_y_continuous(limits = c(min(Genes$logFAP)-1, max(Genes$logFAP)+1)) +
	theme(axis.title=element_text(size=25), axis.text=element_text(size=20), legend.title=element_text(size=25), legend.text=element_text(size=20)) +
	guides(colour= guide_legend(override.aex=list(size=1.5)))
dev.off()

#OnlyRNR
png(filename=paste0(work_dir,"01.Plot/RNR_THBS2_FAP.png"), width=60, height=50, units="cm", res=200)
ggplot(data=OnlyRNR_Genes, mapping=aes(x=logTHBS2,y=logFAP, shape=RNR, color=RNR)) + 
    geom_point(size=5) + 
    geom_smooth(method=lm, se=FALSE) + 
    xlab(expression("log"[2]*"THBS2")) +
    ylab(expression("log"[2]*"FAP")) +
	scale_x_continuous(limits = c(min(Genes$logFAP)-1, max(Genes$logFAP)+1))+
	scale_y_continuous(limits = c(min(Genes$logFAP)-1, max(Genes$logFAP)+1)) +
    theme(axis.title=element_text(size=25), axis.text=element_text(size=20), legend.title=element_text(size=25), legend.text=element_text(size=20)) +
    guides(colour= guide_legend(override.aex=list(size=1.5)))
dev.off()


png(filename=paste0(work_dir,"01.Plot/FAP_IL6.png"), width=60, height=50, units="cm", res=200)
ggplot(Genes, aes(x=logFAP,y=logIL6, shape=response, color=response)) + 
    geom_point(size=5) +
    geom_smooth(method=lm, se=FALSE) +
    xlab(expression("log"[2]*"FAP")) +
    ylab(expression("log"[2]*"IL6")) +
    scale_x_continuous(limits = c(min(Genes$logFAP)-1, max(Genes$logFAP)+1))+
    scale_y_continuous(limits = c(min(Genes$logIL6)-1, max(Genes$logIL6)+1)) +
    theme(axis.title=element_text(size=25), axis.text=element_text(size=20), legend.title=element_text(size=25), legend.text=element_text(size=20)) +
    guides(colour= guide_legend(override.aex=list(size=1.5)))
dev.off()

#OnlyRNR
png(filename=paste0(work_dir,"01.Plot/RNR_FAP_IL6.png"), width=60, height=50, units="cm", res=200)
ggplot(OnlyRNR_Genes, aes(x=logFAP,y=logIL6, shape=RNR, color=RNR)) +
    geom_point(size=5) +
    geom_smooth(method=lm, se=FALSE) +
    xlab(expression("log"[2]*"FAP")) +
    ylab(expression("log"[2]*"IL6")) +
    scale_x_continuous(limits = c(min(Genes$logFAP)-1, max(Genes$logFAP)+1))+
    scale_y_continuous(limits = c(min(Genes$logIL6)-1, max(Genes$logIL6)+1)) +
    theme(axis.title=element_text(size=25), axis.text=element_text(size=20), legend.title=element_text(size=25), legend.text=element_text(size=20)) +
    guides(colour= guide_legend(override.aex=list(size=1.5)))
dev.off()

