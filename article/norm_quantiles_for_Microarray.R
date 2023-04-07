suppressMessages({
library(edgeR)
library(limma)
library(RColorBrewer)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)
library(repr)
library(statmod)
library(GO.db)
library(ComplexHeatmap)
library(mixOmics)
library(dplyr)
library(preprocessCore)
library(hgu133plus2.db)
})

work_dir <- "/data/keeyoung/gy_RNA/GSE16879/"
name="GSE16879"
if (!dir.exists(paste0(work_dir,"01.Plot"))){
    dir.create(paste0(work_dir,"01.Plot"))
}
if (!dir.exists(paste0(work_dir,"02.Table"))){
    dir.create(paste0(work_dir,"02.Table"))
}
if (!dir.exists(paste0(work_dir,"03.Pathway"))){
    dir.create(paste0(work_dir,"03.Pathway"))
}

test <- read.table("/data/keeyoung/gy_RNA/GSE16879/GSE16879_series_matrix_nohead.txt", sep='\t', row.names=1)
colnames(test) <- test[1,]
test <- test[-1,]
#as.numeric
for ( A in 1:length(colnames(test))){
test[,A] <- as.numeric(test[,A])
}
test <- log(test, base=2)

#plot
m <- matrix(ncol=2)

for (A in 1:(length(colnames(test))-1)){

TEST <- as.data.frame(test[,A])
NAME <-rep(colnames(test)[A],54675)
mtx <-cbind(TEST[,1], NAME)
dim(mtx)
m <- rbind(m,mtx)
dim(m)
}

m <- m[-1,]
colnames(m) <- c("expression","Sample")
m <- as.data.frame(m)
m$expression <- as.numeric(m$expression)

pdf("test.pdf")
p <- ggplot(m, aes(x=Sample, y=expression)) + geom_boxplot(outlier.shape=NA)
print({p})
dev.off()

#Quantile norm
df <- test
df_norm <- as.data.frame(normalize.quantiles(as.matrix(df)))
m_norm <- normalize.quantiles(as.matrix(df))
colnames(m_norm) <- colnames(test)
rownames(m_norm) <- rownames(test)

#plot
df_m <- matrix(ncol=2)

for (A in 1:(length(colnames(df_norm))-1)){

df_TEST <- as.data.frame(df_norm[,A])
df_NAME <-rep(colnames(df_norm)[A],54675)
df_mtx <-cbind(df_TEST[,1], df_NAME)
dim(df_mtx)  
df_m <- rbind(df_m,df_mtx)
dim(df_m)
}

df_m <- df_m[-1,]
colnames(df_m) <- c("expression","Sample")
df_m <- as.data.frame(df_m)
df_m$expression <- as.numeric(df_m$expression)

pdf("df_norm.pdf")
q <- ggplot(df_m, aes(x=Sample, y=expression)) + geom_boxplot(outlier.shape=NA)
print({q})
dev.off()

#filter
mode(m_norm) <- "numeric"
gSymbols <- AnnotationDbi::select(hgu133plus2.db, keys=rownames(m_norm), columns="SYMBOL", keytype="PROBEID")
gSymbolsUnique <- gSymbols[!is.na(gSymbols[,2]),]
gSymbolsUnique <- gSymbolsUnique[table(gSymbolsUnique[,1])==1,]
gSymbolsUnique <- sort(setNames(gSymbolsUnique[,2], gSymbolsUnique[,1]))
m_norm <- m_norm[names(gSymbolsUnique),]
rownames(m_norm) <- unname(gSymbolsUnique)

#mean
norm_mtx <- data.frame(matrix(nrow=21855))
for (A in 1:length(colnames(m_norm))){
print(A)
norm_mtx[,A] <- (tapply(m_norm[,A], factor(rownames(m_norm)), mean))
}
colnames(norm_mtx) <- colnames(m_norm)
rownames(norm_mtx) <- rownames(tapply(m_norm[,A], factor(rownames(m_norm)), mean))

m_norm <- norm_mtx

#Only Protein Coding Genes
print("Loading data")
gtf <- read.table("/data/keeyoung/gy_RNA/REFERENCE/gtf/anno_Homo_sapiens.GRCh38.100.gtf", sep=" ", head=T)
m <- match(rownames(m_norm), gtf$gene_symbol)
Func <- gtf$function.[m]
m_norm <- cbind(m_norm, Func)
g <- grep("protein_coding", m_norm[,134])
m_norm <- m_norm[g,]

remove <- grep("^LOC", rownames(m_norm), value=T)
m_norm <- m_norm[which(!rownames(m_norm) %in% remove),]
m_norm <- m_norm[,-134]
colSums(is.infinite(as.matrix(m_norm)))

rm_norm <- m_norm
for ( A in 1:length(colnames(rm_norm))){
w <- which(rm_norm[,A] %in% "-Inf")
if(length(w) > 0){
	rm_norm <- rm_norm[-w,]}
}

nonlog_m_norm <- rm_norm^2
write.table(nonlog_m_norm ,"/data/keeyoung/gy_RNA/GSE16879/GSE16879_QN_nonlog_rmdup_rmInf_rmLOC_OnlyProteinCoding.txt",sep='\t', quote=F) #230117

#after
sample.meta <- read.csv("/data/keeyoung/gy_RNA/01.Samplelist/GSE16879_sample.csv", header=F)
v2 <- grep("after",sample.meta$V2)
sample.meta <- sample.meta[v2,]
sample.ID <- sample.meta$V1
sample.pheno <- sample.meta$V2 #DR, NDR
names(sample.pheno) <- sample.ID
m <- match(names(sample.pheno),colnames(nonlog_m_norm))
iCD <- nonlog_m_norm[,m]

if(all.equal(names(sample.pheno),colnames(iCD))=="TRUE"){
sample.group <- factor(sample.pheno)
sample.design <- model.matrix(~0 + sample.group)
colnames(sample.design) <- levels(sample.group)
}

my.contrasts <- makeContrasts(class=NDR_after-DR_after, levels=sample.design)
cnt <- iCD
y <- DGEList(cnt)
y <- estimateDisp(y,sample.design)
Ffit <- glmQLFit(y, sample.design)
Ftest <- glmQLFTest(Ffit, contrast=my.contrasts[,"class"])

#DecideTests and Plot
dt <- (decideTests(Ftest, p.value=0.05, adjust.method="fdr", lfc=1))
w <- which(dt==0)
qcFtest <- Ftest[-w,]
write.table(qcFtest$table,paste0(work_dir,"02.Table/",name,"_after_qcFtest_table_David.txt"), quote=F, sep=',')

#before
sample.meta <- read.csv("/data/keeyoung/gy_RNA/01.Samplelist/GSE16879_sample.csv", header=F)
v2 <- grep("before",sample.meta$V2)
sample.meta <- sample.meta[v2,]
sample.ID <- sample.meta$V1
sample.pheno <- sample.meta$V2 #DR, NDR
names(sample.pheno) <- sample.ID
m <- match(names(sample.pheno),colnames(nonlog_m_norm))
iCD <- nonlog_m_norm[,m]

if(all.equal(names(sample.pheno),colnames(iCD))=="TRUE"){
sample.group <- factor(sample.pheno)
sample.design <- model.matrix(~0 + sample.group)
colnames(sample.design) <- levels(sample.group)
}

my.contrasts <- makeContrasts(class=NDR_before-DR_before, levels=sample.design)
cnt <- iCD
y <- DGEList(cnt)
y <- estimateDisp(y,sample.design)
Ffit <- glmQLFit(y, sample.design)
Ftest <- glmQLFTest(Ffit, contrast=my.contrasts[,"class"])

#DecideTests and Plot
dt <- (decideTests(Ftest, p.value=0.05, adjust.method="fdr", lfc=1))
w <- which(dt==0)
qcFtest <- Ftest[-w,]

Ftest$table[which(rownames(Ftest$table) %in% "FAP"),]
write.table(qcFtest$table,paste0(work_dir,"02.Table/",name,"_before_qcFtest_table_David.txt"), quote=F, sep=',')

#After vs before in NDR
sample.meta <- read.csv("/data/keeyoung/gy_RNA/01.Samplelist/GSE16879_sample.csv", header=F)
v2 <- grep("^NDR",sample.meta$V2)
sample.meta <- sample.meta[v2,]
sample.ID <- sample.meta$V1
sample.pheno <- sample.meta$V2 #after before
names(sample.pheno) <- sample.ID
m <- match(names(sample.pheno),colnames(nonlog_m_norm))
iCD <- nonlog_m_norm[,m]

if(all.equal(names(sample.pheno),colnames(iCD))=="TRUE"){
sample.group <- factor(sample.pheno)
sample.design <- model.matrix(~0 + sample.group)
colnames(sample.design) <- levels(sample.group)
}

my.contrasts <- makeContrasts(class=NDR_after-NDR_before, levels=sample.design)
cnt <- iCD
y <- DGEList(cnt)
y <- estimateDisp(y,sample.design)
Ffit <- glmQLFit(y, sample.design)
Ftest <- glmQLFTest(Ffit, contrast=my.contrasts[,"class"])

#DecideTests and Plot
dt <- (decideTests(Ftest, p.value=0.05, adjust.method="fdr", lfc=1))
w <- which(dt==0)
qcFtest <- Ftest[-w,]
write.table(qcFtest$table,paste0(work_dir,"02.Table/",name,"_AvsB_NDR_qcFtest_table_David.txt"), quote=F, sep=',')

#After vs before in DR
sample.meta <- read.csv("/data/keeyoung/gy_RNA/01.Samplelist/GSE16879_sample.csv", header=F)
v2 <- grep("^DR",sample.meta$V2)
sample.meta <- sample.meta[v2,]
sample.ID <- sample.meta$V1
sample.pheno <- sample.meta$V2 #after before
names(sample.pheno) <- sample.ID
m <- match(names(sample.pheno),colnames(nonlog_m_norm))
iCD <- nonlog_m_norm[,m]

if(all.equal(names(sample.pheno),colnames(iCD))=="TRUE"){
sample.group <- factor(sample.pheno)
sample.design <- model.matrix(~0 + sample.group)
colnames(sample.design) <- levels(sample.group)
}

my.contrasts <- makeContrasts(class=DR_after-DR_before, levels=sample.design)
cnt <- iCD  
y <- DGEList(cnt)
y <- estimateDisp(y,sample.design)
Ffit <- glmQLFit(y, sample.design)
Ftest <- glmQLFTest(Ffit, contrast=my.contrasts[,"class"])

#DecideTests and Plot
dt <- (decideTests(Ftest, p.value=0.05, adjust.method="fdr", lfc=1))
w <- which(dt==0)
qcFtest <- Ftest[-w,]
write.table(qcFtest$table,paste0(work_dir,"02.Table/",name,"_AvsB_DR_qcFtest_table_David.txt"), quote=F, sep=',')


