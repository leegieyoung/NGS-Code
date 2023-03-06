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

#Total_CDno <- read.table("/data/keeyoung/gy_RNA/06.output/CD_no/02.Table/CD_no_lcpm_Norm.txt",sep='\t')
#Tru_CDno <- read.table("/data/keeyoung/gy_RNA/06.output/CD_no_Tru/02.Table/CD_no_Tru_lcpm_Norm.txt",sep='\t')

#Tru_CDno <- Tru_CDno[na.omit(match(rownames(Total_CDno), rownames(Tru_CDno))),]
#Total_CDno <- Total_CDno[na.omit(match(rownames(Tru_CDno), rownames(Total_CDno))),]

#Tru_CDno <- Tru_CDno[order(rownames(Tru_CDno)),]
#Total_CDno <- Total_CDno[order(rownames(Total_CDno)),]
#Common <- rownames(Tru_CDno)

#########
print("Make_DIR")
pre_dir <- "/data/keeyoung/gy_RNA/06.output/"
name="CD_in"
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
#Extract
v2 <- grep("CD_infla",sample.meta$V2)
sample.meta <- sample.meta[v2,]
#v4 <- grep("Totalseq", sample.meta$V4)
#sample.meta <- sample.meta[v4,]

sample.ID <- sample.meta$V1
sample.ID <- gsub("-",".", sample.ID)
sample.pheno <- sample.meta$V5
sample.pheno <- gsub("^NDR$","NR",sample.pheno)
sample.pheno <- gsub("^DR$","R",sample.pheno)
#sample.pheno <- substr(sample.pheno,1,5)
names(sample.pheno) <- sample.ID

Unknown <- grep("Unknown",sample.pheno)
if(length(Unknown)!=0){
sample.pheno <- sample.pheno[-Unknown]
sample.ID <- names(sample.pheno)}

if(all.equal(sample.ID, names(sample.pheno))!="TRUE"){
print("")
print("ERROR : Not Equal sample.ID and sample.pheno!!!!!!!!!!!!!!")
print("")}


#sample.group <- factor(sample.pheno, levels=c("NSNP","stricture","penetrating"), order = T)
#sample.design <- model.matrix(~0 + sample.group)
#colnames(sample.design) <- levels(sample.group)

#Sample data
Total_CDin <- read.table("/data/keeyoung/gy_RNA/06.output/CD_in/02.Table/CD_in_lcpm_Norm.txt", sep='\t')
Tru_CDin <- read.table("/data/keeyoung/gy_RNA/06.output/CD_in_Tru/02.Table/CD_in_Tru_lcpm_Norm.txt", sep='\t')

#match
Tru_CDin <- Tru_CDin[na.omit(match(rownames(Total_CDin), rownames(Tru_CDin))),]
Total_CDin <- Total_CDin[na.omit(match(rownames(Tru_CDin), rownames(Total_CDin))),]
all.equal(rownames(Tru_CDin), rownames(Total_CDin))


merge <- cbind(Total_CDin,Tru_CDin)

merge <- merge[,match(names(sample.pheno), colnames(merge))]


#It is working when all.equal(sample.ID,colnames(merge)) is "TRUE".
if(all.equal(names(sample.pheno),colnames(merge))=="TRUE"){
sample.group <- factor(sample.pheno)
sample.design <- model.matrix(~0 + sample.group)
colnames(sample.design) <- levels(sample.group)
}

###################
merge.pca <- pca(t(merge))
pca.variates <- merge.pca$variates[[1]]
pca.variates <- pca.variates[,1:2]
pca.variates <- cbind(pca.variates, sample.pheno)
pca.variates <- as.data.frame(pca.variates)
pca.variates$PC1 <- as.numeric(pca.variates$PC1)
pca.variates$PC2 <- as.numeric(pca.variates$PC2)
p <- ggplot(pca.variates, aes(PC1,PC2, shape=factor(sample.pheno)))
p <- p + geom_point(aes(colour = factor(sample.pheno)), size=5)
p <- p + scale_color_manual(values=c('#56B4E9','#E69F00'))

png(filename=paste0(work_dir,"01.Plot/test_pca.png"), width=20, height=15, units="cm", res=200)
print({p})
dev.off()
write.table(pca.variates, paste0(work_dir,"02.Table/test_pca.txt"), sep='\t', quote=F)

train.X <- t(merge)
train.Y <- sample.group
#TRAIN sPLSDA
train.splsda <- splsda(train.X, train.Y, ncomp = 10)
perf.splsda.train <- perf(train.splsda, validation = "Mfold", 
                          folds = 10, nrepeat = 50, # use repeated cross-validation
                          progressBar = TRUE, auc = TRUE, # include AUC values
#						  cpus=36)
						  )

png(filename=paste0(work_dir,"01.Plot/sPLSDA_perf.png"), width=60, height=50, units="cm", res=200)
print({plot(perf.splsda.train, col = color.mixo(5:7), sd = TRUE,
    legend.position = "horizontal")})
dev.off()
saveRDS(perf.splsda.train,paste0(work_dir,"perf.splsda.train.rds"))

perf.splsda.train$choice.ncomp

list.keepX <- c(1:10,  seq(20, 300, 10))

tune.splsda.train <- tune.splsda(train.X, train.Y, ncomp = 1, # calculate for first 4 components
                                 validation = 'Mfold',
                                 folds = 10, nrepeat = 50, # use repeated cross-validation
                                 dist = 'max.dist', # use max.dist measure
                                 measure = "overall", # use balanced error rate of dist measure
                                 test.keepX = list.keepX,
								 progressBar = TRUE,
                                 cpus = 32) # allow for paralleliation to decrease runtime
#)
png(filename=paste0(work_dir,"01.Plot/sPLSDA_tune.png"), width=60, height=50, units="cm", res=200)
print({plot(tune.splsda.train, col = color.jet(1))})
dev.off()
saveRDS(tune.splsda.train,paste0(work_dir,"tune.splsda.train.rds"))

#optimal.ncomp <- tune.splsda.train$choice.ncomp$ncomp
optimal.ncomp <- 1
#optimal.keepX <- tune.splsda.train$choice.keepX[1:optimal.ncomp]
optimal.keepX <- tune.splsda.train$choice.keepX[1]

opti <- cbind(optimal.ncomp, optimal.keepX)
write.table(opti, paste0(work_dir,"optimal.train_keepX.txt"),sep='\t', quote=F)

train.final.splsda <- splsda(train.X, train.Y,
                        ncomp = optimal.ncomp,
                        keepX = optimal.keepX)

train.col <- train.Y
#train.col <- gsub("NSNP","#388ECC",train.col) #blue
train.col <- gsub("^NR$","#C2C2C2",train.col) #gray
train.col <- gsub("^R$","#F68B33",train.col) #orange

legend=list(legend=levels(train.Y),
			col=c("#C2C2C2","#F68B33"),
#            col=c("#F68B33","#C2C2C2", "#388ECC"),
            title="Anti-TNF",
            cex=1.4)

png(filename=paste0(work_dir,"01.Plot/sPLSDA_train_final_heatmap_comp1.png"), width=70, height=30, units="cm", res=200)
print({cim(train.final.splsda, row.sideColors = train.col,
            keysize=c(0.2,0.9),
            legend = legend)
    })
dev.off()

#Loading genes
zero <- (which(train.final.splsda$loadings$X=="0"))
col1 <- rownames(train.final.splsda$loadings$X)[-zero]
col2 <- as.numeric(train.final.splsda$loadings$X[-zero])
train.gene <- as.data.frame(cbind(col1,col2))
train.gene$col2 <- as.numeric(train.gene$col2)
train.gene <- train.gene[order(train.gene[,2], decreasing=TRUE),]
write.table(train.gene, paste0(work_dir,"02.Table/train.gene.txt"),sep='\t', quote=F)

#predict
#test.final.splsda <- predict(train.final.splsda, test.X, dist="max.dist")
#test.max.dist <- test.final.splsda$class$max.dist[,1]
#table(factor(test.max.dist, levels=levels(sample.group)), test.Y)
