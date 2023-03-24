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

#########
print("Make_DIR")
pre_dir <- "/data/keeyoung/gy_RNA/06.output/"
name="CD_no_2group_platform"
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

#sample clinical
sample.meta <- read.table("/data/keeyoung/scRNA/iCD/01.Samplelist/sample.csv", sep=',')
##Extract
	#CD_normal
v2 <- grep("CD_normal",sample.meta$V2)
sample.meta <- sample.meta[v2,]

sample.ID <- sample.meta$V1
sample.ID <- gsub("-",".", sample.ID)
sample.pheno <- sample.meta$V6
names(sample.pheno) <- sample.ID

	#Last_bahavior
Unknown <- grep("Unknown",sample.pheno)
if(length(Unknown)!=0){
sample.pheno <- sample.pheno[-Unknown]}

sample.pheno <- gsub("stricture","other",sample.pheno)
sample.pheno <- gsub("penetrating","other",sample.pheno)

	#platform
sample.platform <- sample.meta$V4
names(sample.platform) <- sample.meta$V1
names(sample.platform) <- gsub("-",".", names(sample.platform))
sample.platform <- sample.platform[which(names(sample.platform) %in% names(sample.pheno))]


#Common genes
Total_CDno <- read.table("/data/keeyoung/gy_RNA/06.output/CD_no/02.Table/CD_no_Common_lcpm_Norm.txt", sep='\t')
Tru_CDno <- read.table("/data/keeyoung/gy_RNA/06.output/CD_no_Tru/02.Table/CD_no_Tru_Common_lcpm_Norm.txt", sep='\t')

sample.group <- factor(sample.pheno, levels=c("NSNP","other"), order = T)
sample.design <- model.matrix(~0 + sample.group)
colnames(sample.design) <- levels(sample.group)

	#match
Tru_CDno <- Tru_CDno[na.omit(match(rownames(Total_CDno), rownames(Tru_CDno))),]
Total_CDno <- Total_CDno[na.omit(match(rownames(Tru_CDno), rownames(Total_CDno))),]
all.equal(rownames(Tru_CDno), rownames(Total_CDno))
CommonGenes <- rownames(Total_CDno)

#############Raw data
print("Loading raw data")
iCD <- read.table("/data/keeyoung/gy_RNA/02.featureCounts/gene_id/454.IBD.fc_P_M.merge",sep='\t', head=T)
colnames(iCD) <- gsub(".TR","",colnames(iCD))

#######################################################
m <- match(names(sample.pheno),colnames(iCD))
iCD <- iCD[,m]
print({
	paste0("The iCD's Dimension is : ", dim(iCD)[2], ", And The number of Samples Dimension is : ", length(sample.group))
	})

#match
iCD <- iCD[,na.omit(match(names(sample.pheno), colnames(iCD)))]
sample.pheno <- sample.pheno[na.omit(match(names(sample.pheno), colnames(iCD)))]

print({   
    paste0("The iCD's Dimension is : ", dim(iCD)[2], ", And The number of Samples Dimension is : ", length(sample.group))
    })


##################
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
print({
	head(y)
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
write.table(y_uniq_cpm_NoQC, paste0(work_dir,"02.Table/",name,"_NoQC_cpm.txt"), quote=F, sep='\t')

print({
    paste0("Extract Protein Codind Genes and Remove LOC genes : ", dim(y_uniq))
    })
#Normalized
keep <- filterByExpr(y_uniq)
y_uniq <- y_uniq[keep, , keep.lib.sizes=FALSE]
y_uniq <- y_uniq[filterByExpr(y_uniq),]
y_uniq <- calcNormFactors(y_uniq)
y_uniq <- estimateDisp(y_uniq,sample.design)
#y_uniq <- estimateDisp(y_uniq) #Using classic mode
print({
    paste0("Do Norm, and gene length is :" ,dim(y_uniq)[1])
    })

y_uniq.lcpm <- cpm(y_uniq, log=TRUE)

###################
y_uniq.lcpm.t <- t(y_uniq.lcpm)
Total.lcpm.t <- y_uniq.lcpm.t[rownames(y_uniq.lcpm.t) %in% names(sample.platform[which(sample.platform %in% "Totalseq")]),]
Tru.lcpm.t <- y_uniq.lcpm.t[rownames(y_uniq.lcpm.t) %in% names(sample.platform[which(sample.platform %in% "Truseq")]),]

Total.pca <- pca(Total.lcpm.t, ncomp=10, center=TRUE, scale=FALSE)
pca.variates <- Total.pca$variates[[1]]
pca.variates <- pca.variates[,1:2]
pca.variates <- cbind(pca.variates, sample.pheno)
pca.variates <- as.data.frame(pca.variates)
pca.variates$PC1 <- as.numeric(pca.variates$PC1)
pca.variates$PC2 <- as.numeric(pca.variates$PC2)
p <- ggplot(pca.variates, aes(PC1,PC2, shape=factor(sample.pheno)))
p <- p + geom_point(aes(colour = factor(sample.pheno)), size=5)
p <- p + scale_color_manual(values=c('#56B4E9','#E69F00'))

png(filename=paste0(work_dir,"01.Plot/filter_Total_pca.png"), width=20, height=15, units="cm", res=200)
print({p})
dev.off()
write.table(pca.variates, paste0(work_dir,"02.Table/filter_Total_pca.txt"), sep='\t', quote=F)

Tru.pca <- pca(Tru.lcpm.t, ncomp=10, center=TRUE, scale=FALSE)
pca.variates <- Tru.pca$variates[[1]]
pca.variates <- pca.variates[,1:2]
pca.variates <- cbind(pca.variates, sample.pheno)
pca.variates <- as.data.frame(pca.variates)
pca.variates$PC1 <- as.numeric(pca.variates$PC1)
pca.variates$PC2 <- as.numeric(pca.variates$PC2)
p <- ggplot(pca.variates, aes(PC1,PC2, shape=factor(sample.pheno)))
p <- p + geom_point(aes(colour = factor(sample.pheno)), size=5)
p <- p + scale_color_manual(values=c('#56B4E9','#E69F00'))

png(filename=paste0(work_dir,"01.Plot/filter_Tru_pca.png"), width=20, height=15, units="cm", res=200)
print({p})
dev.off()
write.table(pca.variates, paste0(work_dir,"02.Table/filter_Tru_pca.txt"), sep='\t', quote=F)

#######
print("Run glmQLF")
#glmQLF = F-test
my.contrasts <- makeContrasts(class=NSNP-other, levels=sample.design)
Ffit <- glmQLFit(y_uniq, sample.design)
Ftest <- glmQLFTest(Ffit, contrast=my.contrasts[,"class"])
print("Ftest")
print({dim(Ftest)})
png(filename=paste0(work_dir,"01.Plot/glmQLF.png"), width=60, height=50, units="cm", res=200)
print({plotMD(Ftest) + abline(h=c(-1, 1), col="blue")})
dev.off()
write.table(Ftest$table,paste0(work_dir,"02.Table/",name,"_Ftest_table.txt"), quote=F, sep='\t')

#DecideTests and Plot
dt <- (decideTests(Ftest, p.value=0.05, adjust.method="fdr", lfc=1))
w <- which(dt==0)
qcFtest <- Ftest[-w,]
print({
    paste0("qcFtest's demension is : ", dim(qcFtest))
    })
png(filename=paste0(work_dir,"01.Plot/QC_glmQLF_L.png"), width=60, height=50, units="cm", res=200)
print({plotMD(qcFtest) + abline(h=c(-1, 1), col="blue")})
dev.off()
write.table(qcFtest$table,paste0(work_dir,"02.Table/",name,"_qcFtest_table.txt"), quote=F, sep='\t')


w <- which(rownames(y_uniq) %in% rownames(qcFtest$genes))
y_uniq_dt <- y_uniq[w,]
y_uniq_dt.lcpm <- cpm(y_uniq_dt, log=TRUE)

##########PCA_dt##########
y_uniq_dt.lcpm.t <- t(y_uniq_dt.lcpm)
Total.lcpm.t <- y_uniq_dt.lcpm.t[rownames(y_uniq_dt.lcpm.t) %in% names(sample.platform[which(sample.platform %in% "Totalseq")]),]
Tru.lcpm.t <- y_uniq_dt.lcpm.t[rownames(y_uniq_dt.lcpm.t) %in% names(sample.platform[which(sample.platform %in% "Truseq")]),]

Total.pca <- pca(Total.lcpm.t, ncomp=10, center=TRUE, scale=FALSE)
pca.variates <- Total.pca$variates[[1]]
pca.variates <- pca.variates[,1:2]
pca.variates <- cbind(pca.variates, sample.pheno)
pca.variates <- as.data.frame(pca.variates)
pca.variates$PC1 <- as.numeric(pca.variates$PC1)
pca.variates$PC2 <- as.numeric(pca.variates$PC2)
p <- ggplot(pca.variates, aes(PC1,PC2, shape=factor(sample.pheno)))
p <- p + geom_point(aes(colour = factor(sample.pheno)), size=5)
p <- p + scale_color_manual(values=c('#56B4E9','#E69F00'))

png(filename=paste0(work_dir,"01.Plot/dt_Total_pca.png"), width=20, height=15, units="cm", res=200)
print({p})
dev.off()
write.table(pca.variates, paste0(work_dir,"02.Table/dt_Total_pca.txt"), sep='\t', quote=F)

Tru.pca <- pca(Tru.lcpm.t, ncomp=10, center=TRUE, scale=FALSE)
pca.variates <- Tru.pca$variates[[1]]
pca.variates <- pca.variates[,1:2]
pca.variates <- cbind(pca.variates, sample.pheno)
pca.variates <- as.data.frame(pca.variates)
pca.variates$PC1 <- as.numeric(pca.variates$PC1)
pca.variates$PC2 <- as.numeric(pca.variates$PC2)
p <- ggplot(pca.variates, aes(PC1,PC2, shape=factor(sample.pheno)))
p <- p + geom_point(aes(colour = factor(sample.pheno)), size=5)
p <- p + scale_color_manual(values=c('#56B4E9','#E69F00'))

png(filename=paste0(work_dir,"01.Plot/dt_Tru_pca.png"), width=20, height=15, units="cm", res=200)
print({p})
dev.off()
write.table(pca.variates, paste0(work_dir,"02.Table/dt_Tru_pca.txt"), sep='\t', quote=F)



###################
Total.sample <- sample.platform[which(sample.platform %in% "Totalseq")]
Tru.sample <- sample.platform[which(sample.platform %in% "Truseq")]

train.X <- y_uniq_dt.lcpm[,which(colnames(y_uniq_dt.lcpm) %in% names(Total.sample))]
test.X <- y_uniq_dt.lcpm[,which(colnames(y_uniq_dt.lcpm) %in% names(Tru.sample))]

train.Y <- sample.group[which(names(sample.group) %in% names(Total.sample))]
test.Y <- sample.group[which(names(sample.group) %in% names(Tru.sample))]

if(all.equal(colnames(train.X), names(train.Y))=="TRUE"){
	if(all.equal(colnames(test.X), names(test.Y))=="TRUE"){
train.X <- t(train.X)
test.X <- t(test.X)
}}

#TRAIN sPLSDA
train.splsda <- splsda(train.X, train.Y, ncomp = 10)
perf.splsda.train <- perf(train.splsda, validation = "Mfold", 
                          folds = 10, nrepeat = 50, # use repeated cross-validation
                          progressBar = TRUE, auc = TRUE, # include AUC values
						  cpus=36)
#						  )

png(filename=paste0(work_dir,"01.Plot/dt_sPLSDA_perf.png"), width=60, height=50, units="cm", res=200)
print({plot(perf.splsda.train, col = color.mixo(5:7), sd = TRUE,
    legend.position = "horizontal")})
dev.off()
saveRDS(perf.splsda.train,paste0(work_dir,"dt_perf.splsda.train.rds"))
#perf.splsda.train <- readRDS(paste0(work_dir,"dt_perf.splsda.train.rds"))

perf.splsda.train$choice.ncomp
write.table(perf.splsda.train$choice.ncomp, paste0(work_dir,"perf.table.txt"), sep='\t', quote=F)
list.keepX <- c(1:10,  seq(20, 300, 10))

tune.splsda.train <- tune.splsda(train.X, train.Y, ncomp = 2, # calculate for first 4 components
                                 validation = 'Mfold',
                                 folds = 10, nrepeat = 50, # use repeated cross-validation
                                 dist = 'centroids.dist', # use max.dist measure
                                 measure = "overall", # use balanced error rate of dist measure
                                 test.keepX = list.keepX,
								 progressBar = TRUE,
                                 cpus = 32) # allow for paralleliation to decrease runtime
#
png(filename=paste0(work_dir,"01.Plot/dt_sPLSDA_tune.png"), width=60, height=50, units="cm", res=200)
print({plot(tune.splsda.train, col = color.jet(2))})
dev.off()
saveRDS(tune.splsda.train,paste0(work_dir,"dt_tune.splsda.train.rds"))
#tune.splsda.train <- readRDS(paste0(work_dir,"dt_tune.splsda.train.rds"))
optimal.ncomp <- tune.splsda.train$choice.ncomp$ncomp
#optimal.ncomp <- 1
optimal.keepX <- tune.splsda.train$choice.keepX[1:optimal.ncomp]
#optimal.keepX <- tune.splsda.train$choice.keepX[1]

opti <- cbind(optimal.ncomp, optimal.keepX)
write.table(opti, paste0(work_dir,"dt_optimal.train_keepX.txt"),sep='\t', quote=F)

train.final.splsda <- splsda(train.X, train.Y,
                        ncomp = optimal.ncomp,
                        keepX = optimal.keepX)

png(filename=paste0(work_dir,"01.Plot/dt_sPLSDA_pca.png"), width=60, height=50, units="cm", res=200)
print({plotIndiv(train.final.splsda, comp = 1:2,
        group = train.Y, ind.names = FALSE,
        ellipse = TRUE,
        legend = TRUE,
        title="sPLS-DA"
)})
dev.off()


train.col <- train.Y
train.col <- gsub("NSNP","#388ECC",train.col) #blue
train.col <- gsub("other","#F68B33",train.col) #orange

legend=list(legend=c("NSNP","other"),
			col=c("#388ECC","#F68B33"),
#            col=c("#F68B33","#C2C2C2", "#388ECC"),
            title="Last Behavior",
            cex=1.4)

png(filename=paste0(work_dir,"01.Plot/dt_sPLSDA_train_final_heatmap_comp1.png"), width=70, height=30, units="cm", res=200)
print({cim(train.final.splsda, row.sideColors = train.col,
            keysize=c(0.2,0.9),
            legend = legend)
    })
dev.off()

#Loading Plot
png(filename=paste0(work_dir,"01.Plot/dt_sPLSDA_loading_Plot_comp1.png"),width=30, height=30, units="cm", res=200)
print({
plotLoadings(train.final.splsda, comp=1, method = 'median', contrib = 'max', size.name=0.7, ndisplay=20, title=paste0("Loadings n comp 1"))
})
dev.off()


png(filename=paste0(work_dir,"01.Plot/dt_sPLSDA_loading_Plot.png"),width=30, height=30, units="cm", res=200)
print({
plotLoadings(train.final.splsda, comp=optimal.ncomp, method = 'median', contrib = 'max', size.name=0.7, ndisplay=20, title=paste0("Loadings n comp ",optimal.ncomp))
})
dev.off()


#Loading genes
zero <- (which(train.final.splsda$loadings$X=="0"))
col1 <- rownames(train.final.splsda$loadings$X)[-zero]
col2 <- as.numeric(train.final.splsda$loadings$X[-zero])
train.gene <- as.data.frame(cbind(col1,col2))
train.gene$col2 <- as.numeric(train.gene$col2)
train.gene <- train.gene[order(train.gene[,2], decreasing=TRUE),]
write.table(train.gene, paste0(work_dir,"02.Table/dt_train.gene.txt"),sep='\t', quote=F)

#predict
test.final.splsda <- predict(train.final.splsda, test.X, dist="max.dist")
test.max.dist <- test.final.splsda$class$max.dist[,1]
Predict_result <- table(factor(test.max.dist, levels=levels(sample.group)), test.Y)
write.table(Predict_result, paste0(work_dir,"02.Table/dt_Predict_result.txt"),sep='\t', quote=F)

#The Performance Plots
png(filename=paste0(work_dir,"01.Plot/dt_sPLSDA_trian_final_performance_plot.png"),width=20, height=15, units="cm", res=200)
print({
	auc.splsda = auroc(train.final.splsda, roc.comp = optimal.ncomp, print = FALSE)
})
dev.off()

##############
#TEST sPLSDA
test.splsda <- splsda(test.X, test.Y, ncomp = 10)
perf.splsda.test <- perf(test.splsda, validation = "Mfold",
                          folds = 10, nrepeat = 50, # use repeated cross-validation
                          progressBar = TRUE, auc = TRUE, # include AUC values
#                         cpus=36)
                          )

png(filename=paste0(work_dir,"01.Plot/dt_sPLSDA_perf_test.png"), width=60, height=50, units="cm", res=200)
print({plot(perf.splsda.test, col = color.mixo(5:7), sd = TRUE,
    legend.position = "horizontal")})
dev.off()
saveRDS(perf.splsda.test,paste0(work_dir,"dt_perf.splsda.test.rds"))

perf.splsda.test$choice.ncomp

list.keepX <- c(1:10,  seq(20, 300, 10))

tune.splsda.test <- tune.splsda(test.X, test.Y, ncomp = 3, # calculate for first 4 components
                                 validation = 'Mfold',
                                 folds = 10, nrepeat = 50, # use repeated cross-validation
                                 dist = 'max.dist', # use max.dist measure
                                 measure = "overall", # use balanced error rate of dist measure
                                 test.keepX = list.keepX,
                                 progressBar = TRUE,
                                 cpus = 32) # allow for paralleliation to decrease runtime
#)
png(filename=paste0(work_dir,"01.Plot/dt_sPLSDA_tune_test.png"), width=60, height=50, units="cm", res=200)
print({plot(tune.splsda.test, col = color.jet(3))})
dev.off()
saveRDS(tune.splsda.test,paste0(work_dir,"dt_tune.splsda.test.rds"))

#optimal.ncomp <- tune.splsda.test$choice.ncomp$ncomp
optimal.ncomp <- 1
#optimal.keepX <- tune.splsda.test$choice.keepX[1:optimal.ncomp]
optimal.keepX <- tune.splsda.test$choice.keepX[1]

opti <- cbind(optimal.ncomp, optimal.keepX)
write.table(opti, paste0(work_dir,"dt_optimal.test_keepX.txt"),sep='\t', quote=F)

test.final.splsda <- splsda(test.X, test.Y,
                        ncomp = optimal.ncomp,
                        keepX = optimal.keepX)

test.col <- test.Y
test.col <- gsub("NSNP","#388ECC",test.col) #blue
test.col <- gsub("other","#F68B33",test.col) #orange

legend=list(legend=levels(test.Y),
            col=c("#388ECC","#F68B33"),
#            col=c("#F68B33","#C2C2C2", "#388ECC"),
            title="Last Behavior",
            cex=1.4)

png(filename=paste0(work_dir,"01.Plot/dt_sPLSDA_test_final_heatmap_comp1.png"), width=70, height=30, units="cm", res=200)
print({cim(test.final.splsda, row.sideColors = test.col,
            keysize=c(0.2,0.9),
            legend = legend)
    })
dev.off()


#Loading Plot
png(filename=paste0(work_dir,"01.Plot/dt_sPLSDA_test_loading_Plot_comp1.png"),width=30, height=30, units="cm", res=200)
print({
plotLoadings(train.final.splsda, comp=1, method = 'median', contrib = 'max', size.name=0.7, ndisplay=20, title=paste0("Loadings n comp 1"))
})
dev.off()

#Loading genes
zero <- (which(test.final.splsda$loadings$X=="0"))
col1 <- rownames(test.final.splsda$loadings$X)[-zero]
col2 <- as.numeric(test.final.splsda$loadings$X[-zero])
test.gene <- as.data.frame(cbind(col1,col2))
test.gene$col2 <- as.numeric(test.gene$col2)
test.gene <- test.gene[order(test.gene[,2], decreasing=TRUE),]
write.table(test.gene, paste0(work_dir,"02.Table/dr_test.gene.txt"),sep='\t', quote=F)

#test.final.splsda <- predict(test.final.splsda, test.X, dist="max.dist")


