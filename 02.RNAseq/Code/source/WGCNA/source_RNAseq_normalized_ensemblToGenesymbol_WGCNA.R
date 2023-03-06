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
})

RNAseq <- function(name,iCD){
#DIR
source("/data/keeyoung/gy_RNA/Code/source/source_RNAseq_go.R")
source("/data/keeyoung/gy_RNA/Code/source/source_RNAseq_kegg.R")

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
if (!dir.exists(paste0(work_dir,"04.WGCNA"))){
    dir.create(paste0(work_dir,"04.WGCNA"))
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
	paste0("Do Norm is :" ,dim(y_uniq))
	})
y_uniq_cpm <- cpm(y_uniq, log=F)
write.table(y_uniq_cpm, paste0(work_dir,"02.Table/",name,"_cpm.txt"), quote=F, sep='\t')

#####Analysis####
#PCA
y_uniq.lcpm <- cpm(y_uniq, log=TRUE)
y_uniq.lcpm.t <- t(y_uniq.lcpm)
y_uniq.pca <- pca(y_uniq.lcpm.t, ncomp=10, center=TRUE, scale=FALSE)
png(filename=paste0(work_dir,"01.Plot/allgenes_pca.png"), width=20, height=15, units="cm", res=200)
pca.variates <- y_uniq.pca$variates[[1]]
pca.variates <- pca.variates[,1:2]
pca.variates <- cbind(pca.variates, sample.pheno)
pca.variates <- as.data.frame(pca.variates)
pca.variates$PC1 <- as.numeric(pca.variates$PC1)
pca.variates$PC2 <- as.numeric(pca.variates$PC2)
p <- ggplot(pca.variates, aes(PC1,PC2, shape=factor(sample.pheno)))
p <- p + geom_point(aes(colour = factor(sample.pheno)), size=5)
p <- p + scale_color_manual(values=c('#56B4E9','#E69F00'))
#p < p + scale_x_continuous(limits = c(-150,150)) + scale_y_continuous(limits = c(-150,150))
p <- p + coord_cartesian(xlim = c(min(pca.variates$PC1)-5, max(pca.variates$PC1)+5))
print({p})
dev.off()

print("Run glmQy_uniqLF")
#glmQLF = F-test
#my.contrasts <- makeContrasts(class=NR-R, levels=sample.design)
Ffit <- glmQLFit(y_uniq, sample.design)
Ftest <- glmQLFTest(Ffit, contrast=my.contrasts[,"class"])
print("Ftest")
print({dim(Ftest)})
png(filename=paste0(work_dir,"01.Plot/glmQLF.png"), width=60, height=50, units="cm", res=200)
print({plotMD(Ftest) + abline(h=c(-1, 1), col="blue")})
dev.off()
write.table(Ftest$table,paste0(work_dir,"02.Table/",name,"_Ftest_table.txt"), quote=F, sep='\t')

#Pathway analysis of Ftest genes
print("Pathway analyzed Ftest")
#GO_func(work_dir,Ftest)
#KEGG_func(work_dir,Ftest)

idfound <- Ftest$gene$gene_name %in% mappedRkeys(org.Hs.egSYMBOL)
Ftest <- Ftest[idfound,]
egENTREZID <- toTable(org.Hs.egSYMBOL)
entrez_m <- match(Ftest$genes$gene_name , egENTREZID$symbol)
Ftest$genes$entrezID <- egENTREZID$gene_id[entrez_m]

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

idfound <- qcFtest$gene$gene_name %in% mappedRkeys(org.Hs.egSYMBOL)
qcFtest <- qcFtest[idfound,]
egENTREZID <- toTable(org.Hs.egSYMBOL)
entrez_m <- match(qcFtest$genes$gene_name , egENTREZID$symbol)
qcFtest$genes$entrezID <- egENTREZID$gene_id[entrez_m]
David <- cbind(qcFtest$genes$entrezID ,qcFtest$table$logFC, qcFtest$gene$gene_name)
print({
    paste0("David's dimention is : ", dim(David))
    })
write.table(David,paste0(work_dir,"02.Table/",name,"_qcFtest_table_David.txt"), quote=F, sep=',')

############WGCNA#####################
#00.input data
datExpr <- t(cpm(y_uniq[which(rownames(y_uniq) %in% rownames(qcFtest)),], log=T))

#Select threshold through Scale-free Network created with expression data.
powers=c(c(1:50))
#sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose=5)

pdf(paste0(work_dir,"04.WGCNA/Scale independence.pdf"), width=10, height=6)
cex1 = 0.9;
print({plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"))})
print({text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red")})
print({abline(h=0.90,col="red")})
dev.off()

pdf(paste0(work_dir,"04.WGCNA/Mean connectivity.pdf"), width=10, height=6)
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

softPower=23
adjacency = adjacency(datExpr, power = softPower)
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method = "average");

minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
	deepSplit = 2, pamRespectsDendro = FALSE,
	minClusterSize = minModuleSize);

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
pdf(paste0(work_dir,"04.WGCNA/Pre_dengrogram.pdf"), width=10, height=12)
print({plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
	dendroLabels = FALSE, hang = 0.03,
	addGuide = TRUE, guideHang = 0.05,
	main = "Gene dendrogram and module colors")})
dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf(paste0(work_dir,"04.WGCNA/Cluster_module_eigengenes.pdf"), width=5, height=6)
print({plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")})
th1=0.25
th2=0.15
th3=0.1
# Plot the cut line into the dendrogram
abline(h=th1, col = "red")
abline(h=th2, col = "blue")
abline(h=th3, col = "green")
dev.off()

th1_merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = th1, verbose = 3)
# The merged module colors
th1_mergedColors = th1_merge$colors;
# Eigengenes of the new merged modules:
th1_mergedMEs = th1_merge$newMEs;

th2_merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = th2, verbose = 3)
# The merged module colors
th2_mergedColors = th2_merge$colors;
# Eigengenes of the new merged modules:
th2_mergedMEs = th2_merge$newMEs;

th3_merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = th3, verbose = 3)
# The merged module colors
th3_mergedColors = th3_merge$colors;
# Eigengenes of the new merged modules:
th3_mergedMEs = th3_merge$newMEs;

length(unique(sort(th1_mergedColors)))
length(unique(sort(th2_mergedColors)))
length(unique(sort(th3_mergedColors)))

pdf(paste0(work_dir,"04.WGCNA/After_dengrogram.pdf"), width=10, height=12)
plotDendroAndColors(geneTree, cbind(dynamicColors, th1_mergedColors, th2_mergedColors, th3_mergedColors),
	c("Dynamic Tree Cut", "th1","th2","th3"),
	dendroLabels = FALSE, hang = 0.03,
	addGuide = TRUE, guideHang = 0.05)
dev.off()

MEDissThres = th3

merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


moduleColors = mergedColors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

datTraits <- as.data.frame(sample.pheno)
datTraits$sample.pheno <- gsub("L_NR",2,datTraits$sample.pheno)
datTraits$sample.pheno <- gsub("L_R",1,datTraits$sample.pheno)
datTraits$sample.pheno <- as.numeric(datTraits$sample.pheno)

moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

pdf(paste0(work_dir,"04.WGCNA/Module-trait_relationships.pdf"), width=10, height=6)
labeledHeatmap(Matrix = moduleTraitCor,
	xLabels = names(datTraits),
	yLabels = names(MEs),
	ySymbols = names(MEs),
	colorLabels = FALSE,
	colors = greenWhiteRed(50),
	textMatrix = textMatrix,
	setStdMargins = FALSE,
	cex.text = 0.5,
	zlim = c(-1,1),
	main = paste("Module-trait relationships"))
dev.off()

Cell_frac <- read.table("/data/keeyoung/gy_RNA/01.Samplelist/MGI_RNR_CD_in.txt", sep='\t', head=T)
rownames(Cell_frac) <- Cell_frac[,1]
Cell_frac <- Cell_frac[,-1]

moduleTraitCor = cor(MEs, Cell_frac, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

pdf(paste0(work_dir,"04.WGCNA/Cell_frac_Module-trait_relationships.pdf"), width=10, height=6)
labeledHeatmap(Matrix = moduleTraitCor,
    xLabels = names(Cell_frac),
    yLabels = names(MEs),
    ySymbols = names(MEs),
    colorLabels = FALSE,
    colors = greenWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 0.5,
    zlim = c(-1,1),
    main = paste("Module-trait relationships"))
dev.off()

# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = softPower);
plotTOM = dissTOM^7;
diag(plotTOM) = NA;

pdf(paste0(work_dir,"04.WGCNA/Network_heatmap_plot_all_genes.pdf"), width=10, height=8)
print({TOMplot(plotTOM, geneTree, moduleColors, terrainColors=TRUE,  main = "Network heatmap plot, all genes")})
dev.off()

MET = orderMEs(cbind(MEs, L_DRNDR))

pdf(paste0(work_dir,"04.WGCNA/Eigengene_adjacency_heatmap.pdf"), width=10, height=8)
print({plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)})
dev.off()

genes <- colnames(datExpr)
names(genes) <- moduleColors
genes[which(genes %in% c("FAP","THBS2","TNF","STAT4"))]

#Chapter3

GS1=as.numeric(cor(y.pheno,datExpr, use="p"))
GeneSignificance=abs(GS1)
# Next module significance is defined as average gene significance.
ModuleSignificance=tapply(GeneSignificance, moduleColors, mean, na.rm=T)


}
