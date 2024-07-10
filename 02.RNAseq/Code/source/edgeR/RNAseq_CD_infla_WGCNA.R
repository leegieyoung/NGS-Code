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
library(WGCNA)
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
#sample.pheno <- sample.meta$V2 #Only_CD_in
sample.pheno <- sample.meta$V3 #Likely DR,NDR
sample.pheno <- gsub("^L$","L_R",sample.pheno)
sample.pheno <- gsub("^R$","L_NR",sample.pheno)

sample.pheno <- substr(sample.pheno,1,5)
names(sample.pheno) <- sample.ID

sample.group <- factor(sample.pheno)

m <- match(names(sample.pheno),colnames(iCD))
iCD <- iCD[,m]

if(all.equal(names(sample.pheno),colnames(iCD))=="TRUE"){
sample.group <- factor(sample.pheno)
sample.design <- model.matrix(~0 + sample.group)
colnames(sample.design) <- levels(sample.group)
}

my.contrasts <- makeContrasts(class=L_NR-L_R, levels=sample.design)

#Analysis
name="CD_in_Total_article"
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

symbol_gtf$Length <- as.numeric(symbol_gtf$Length)
symbol_gtf$Start <- as.numeric(symbol_gtf$Start)
symbol_gtf$End <- as.numeric(symbol_gtf$End)

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

y_uniq.lcpm <- cpm(y_uniq, log=TRUE)

print("Run glmQy_uniqLF")
#glmQLF = F-test
#my.contrasts <- makeContrasts(class=NR-R, levels=sample.design)
Ffit <- glmQLFit(y_uniq, sample.design)
Ftest <- glmQLFTest(Ffit, contrast=my.contrasts[,"class"])

idfound <- Ftest$gene$gene_name %in% mappedRkeys(org.Hs.egSYMBOL)
Ftest <- Ftest[idfound,]
egENTREZID <- toTable(org.Hs.egSYMBOL)
entrez_m <- match(Ftest$genes$gene_name , egENTREZID$symbol)
Ftest$genes$entrezID <- egENTREZID$gene_id[entrez_m]

#DecideTests and Plot
dt <- (decideTests(Ftest, p.value=0.05, adjust.method="fdr", lfc=1))
w <- which(dt==0)
qcFtest <- Ftest[-w,]
idfound <- qcFtest$gene$gene_name %in% mappedRkeys(org.Hs.egSYMBOL)
qcFtest <- qcFtest[idfound,]

#WGCNA R vs NR
#The qcFtest
#y_uniq.lcpm <-y_uniq.lcpm[which(rownames(y_uniq.lcpm) %in% rownames(qcFtest)),]

#The 246 genes
y_uniq.lcpm <- cpm(y_uniq, log=TRUE)
valified_genes <- scan("/data/keeyoung/gy_RNA/Code/source/246_gene.list", what=character(0))
valified_genes <- c(valified_genes,"TNF","STAT4")
y_uniq.lcpm <-y_uniq.lcpm[which(rownames(y_uniq.lcpm) %in% valified_genes),]

#The up-regulated
#Upregulated <- rownames(qcFtest[qcFtest$table$logFC > 1,])
#y_uniq.lcpm <- cpm(y_uniq, log=TRUE)
#y_uniq.lcpm <-y_uniq.lcpm[which(rownames(y_uniq.lcpm) %in% Upregulated),]

#QC
datExpr <-t(y_uniq.lcpm[,which(colnames(y_uniq.lcpm) %in% names(sample.pheno[sample.pheno %in% c("L_NR","L_R")]))])
y <- sample.pheno[sample.pheno %in% c("L_NR","L_R")]
y <- gsub("L_NR",2,y)
y <- gsub("L_R",1,y)
y <- as.numeric(y)

meanExpressionByArray=apply( datExpr,1,mean, na.rm=T)
NumberMissingByArray=apply( is.na(data.frame(datExpr)),1, sum)

KeepArray= NumberMissingByArray<500
table(KeepArray)
datExpr=datExpr[KeepArray,]
y=y[KeepArray]

#ADJ1
ADJ1=abs(cor(datExpr,use="p"))^6
k=as.vector(apply(ADJ1,2,sum, na.rm=T)) # When you have relatively few genes (<5000) use the following code
#k=softConnectivity(datE=datExpr,power=6) # When you have a lot of genes use the following code
datExpr=datExpr[, rank(-k,ties.method="first" )<=3600] #For computational reasons, we restrict the network analysis to 3600 most connected genes.
ADJ1=abs(cor(datExpr,use="p"))^6

#TOM(Topological Overlap Matrix) and color clustering
dissTOM=TOMdist(ADJ1)
hierTOM = hclust(as.dist(dissTOM),method="average")
#We now use two Dynamic Tree Cut methods in which the height cut-off plays a minor role. The first method is 
#called the “tree” method and only uses the dendrogram as input.

#The second method is called “hybrid” and is a hybrid between hclust and pam. As input it requires both a dendrogram
#and the dissimilarity that was used to create the dendrogram.
colorDynamicHybridTOM = labels2colors(cutreeDynamic(hierTOM, distM= dissTOM , cutHeight = 0.998,
    deepSplit=2, pamRespectsDendro = FALSE))
colorh1= colorDynamicHybridTOM

#datME
datME=moduleEigengenes(datExpr,colorh1)$eigengenes
colorME=moduleEigengenes(datExpr,colorh1)$validColors

#Cluster_List
CL <- vector(mode="list",length=length(unique(sort(colorh1))))
for (A in 1:length(unique(sort(colorh1)))){
	CL[[A]] <- colnames(datExpr)[which(colorh1 %in% unique(sort(colorh1))[A])]
}

#NS1, Must as numeric of y, like DEGs
NS1=networkScreening(y=y, datME=datME, datExpr=datExpr,
    oddPower=3, blockSize=1000, minimumSampleSize=4,
    addMEy=TRUE, removeDiag=FALSE, weightESy=0.5)

#Top genes
topList=rank(NS1$p.Weighted,ties.method="first")<=30
gene.names= colnames(datExpr)[topList]
#gene.names=c("STAT4","IL6","TNF","THBS2","CXCL5","FAP")
#gene.names=colnames(datExpr[,which(colnames(datExpr) %in% valified_genes)])

#png(filename=paste0(work_dir,"01.Plot/signed correlations.png"), width=60, height=50, units="cm", res=200)
pdf(paste0(work_dir,"01.Plot/signed correlations.pdf"), width=25, height=30)
plotNetworkHeatmap(datExpr, plotGenes = gene.names,
    networkType="signed", useTOM=FALSE,
    power=1, main="signed correlations")
dev.off()

pdf(paste0(work_dir,"01.Plot/unsigned correlations.pdf"), width=25, height=30)
plotNetworkHeatmap(datExpr, plotGenes = gene.names,
    networkType="unsigned", useTOM=FALSE,
    power=1, main="unsigned correlations")
dev.off()

pdf(paste0(work_dir,"01.Plot/signed TOM in a correlations.pdf"), width=25, height=30)
plotNetworkHeatmap(datExpr, plotGenes = gene.names,
    networkType="signed", useTOM=TRUE,
    power=12, main="TOM in a signed network")
dev.off()

pdf(paste0(work_dir,"01.Plot/unsigned TOM in a correlations.pdf"), width=25, height=30)
plotNetworkHeatmap(datExpr, plotGenes = gene.names,
    networkType="unsigned", useTOM=TRUE,
    power=6, main="TOM in an unsigned network")
dev.off()


