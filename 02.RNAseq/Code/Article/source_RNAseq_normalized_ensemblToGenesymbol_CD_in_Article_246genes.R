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
GO_func(work_dir,Ftest)
KEGG_func(work_dir,Ftest)

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

#246 genes table
valified_genes <- scan("/data/keeyoung/gy_RNA/Code/source/246_gene.list", what=character(0))
valified_Ftest <- Ftest[which(rownames(Ftest) %in% valified_genes),]
print({
    paste0("valified_Ftest's demension is : ", dim(valified_Ftest))
	})
png(filename=paste0(work_dir,"01.Plot/glmQLF_L_246.png"), width=30, height=25, units="cm", res=200)
	print({plotMD(valified_Ftest) + abline(h=c(-1, 1), col="blue")})
dev.off()

valified_Ftest_246genes <-topTags(valified_Ftest, n=rownames(valified_Ftest))
write.table(valified_Ftest_246genes,paste0(work_dir,"02.Table/",name,"_valified_Ftest_246genes.txt"), quote=F, sep=',')

#length(rownames(Ftest))
#Ftest
valified_Ftest <-topTags(Ftest, n=length(rownames(Ftest)))
valified_Ftest <- valified_Ftest$table[,c("logFC","FDR")]

#Expression
valified_Ftest <- valified_Ftest %>%
	mutate(
	Expression=case_when(logFC >= log(2) & FDR <= 0.05 ~ "Up-regulated",
  						logFC <= -log(2) & FDR <= 0.05 ~ "Down-regulated",
						TRUE ~ "Unchanged")
)
#anno_246 genes
for(i in 1:length(rownames(valified_Ftest))){
	if(valified_Ftest[,"Expression"][i]=="Up-regulated"){
	if(is.na(match(rownames(valified_Ftest)[i], valified_genes))!="TRUE"){
		valified_Ftest[,"Expression"][i] <- "uniq"}
#	else {
#		valified_Ftest[,"Expression"][i] <- "common"}
	}
}

png(filename=paste0(work_dir,"01.Plot/246gene_volcanoplot.png"), width=22, height=20, units="cm", res=200)
print({ggplot(valified_Ftest, aes(logFC, -log(FDR,10))) +
    geom_point(aes(color =Expression), size = 2/5) +    
    xlab(expression("log"[2]*"FC")) +   
    ylab(expression("-log"[10]*"FDR")) +
    scale_color_manual(values = c("dodgerblue3", "gray50","green","firebrick3")) +
    guides(colour = guide_legend(override.aes = list(size=1.5)))
})
dev.off()


valified_Ftest[,"uniq_genes"] <- ""
for(i in 1:length(rownames(valified_Ftest))){
valified_Ftest[,"uniq_genes"][i] <- valified_genes[match(rownames(valified_Ftest)[i], valified_genes)]
}
valified_Ftest$uniq_genes[is.na(valified_Ftest$uniq_genes)] <- "common"
valified_Ftest$uniq_genes[which(valified_Ftest$uniq_genes!="common")] <- "uniq"

png(filename=paste0(work_dir,"01.Plot/246gene_volcanoplot.png"), width=22, height=20, units="cm", res=200)
print({ggplot(valified_Ftest, aes(logFC, -log(FDR,10))) +
    geom_point(aes(color =uniq_genes), size = 2/5) +    
    xlab(expression("log"[2]*"FC")) +   
    ylab(expression("-log"[10]*"FDR")) +
    scale_color_manual(values = c("dodgerblue3", "gray50","firebrick3")) +
    guides(colour = guide_legend(override.aes = list(size=1.5)))
})
dev.off()

write.table(David,paste0(work_dir,"02.Table/",name,"_qcFtest_table_David.txt"), quote=F, sep=',')

#Genes that pass through fdr and logFC
w <- which(rownames(y_uniq) %in% rownames(qcFtest$genes))
y_uniq_dt <- y_uniq[w,]
y_uniq_dt_heatmap <- y_uniq_dt
rownames(y_uniq_dt_heatmap) <- y_uniq_dt$genes$gene_name
#y_uniq_dt.lcpm <- cpm(y_uniq_dt, log=T)
y_uniq_dt.cpm <- cpm(y_uniq_dt, log=F)
write.table(y_uniq_dt.cpm, paste0(work_dir,"02.Table/",name,"_qcFtest_cpm.txt"), quote=F, sep='\t')

y_uniq_dt.lcpm <- cpm(y_uniq_dt_heatmap, log=T)
print({
	paste0("PCA gene's dimension is ", dim(y_uniq_dt.lcpm))
	})

###heatmap############
#246 Genes 
w_vali <- which(y_uniq$genes$gene_name %in% valified_genes)
y_uniq_vali <- y_uniq[w_vali,]
y_uniq_vali_heatmap <- y_uniq_vali
rownames(y_uniq_vali_heatmap) <- y_uniq_vali$genes$gene_name
y_uniq_vali.lcpm <- cpm(y_uniq_vali_heatmap, log=T)

#IBD_inflammatory_fibroblasts, Iif
Iif_genes <- c("STAT4","TNF","IL17","IRF1","IL6","FAP")
wIif <- which(y_uniq$genes$gene_name %in% Iif_genes)
y_uniq_Iif <- y_uniq[wIif,]
y_uniq_Iif_heatmap <- y_uniq_Iif
rownames(y_uniq_Iif_heatmap) <- y_uniq_Iif$genes$gene_name
y_uniq_Iif.lcpm <- cpm(y_uniq_Iif_heatmap, log=T)

test <- as.data.frame(t(y_uniq_Iif.lcpm))

TNF_IL6_R <- cor.test(x=test[which(rownames(test) %in% names(sample.pheno[which(sample.pheno %in% "L_R")])),][,"TNF"],y=test[which(rownames(test) %in% names(sample.pheno[which(sample.pheno %in% "L_R")])),][,"IL6"],method="pearson")
print({"TNF_IL6_R"})
print({TNF_IL6_R})
TNF_IL6_NR <- cor.test(x=test[which(rownames(test) %in% names(sample.pheno[which(sample.pheno %in% "L_NR")])),][,"TNF"],y=test[which(rownames(test) %in% names(sample.pheno[which(sample.pheno %in% "L_NR")])),][,"IL6"],method="pearson")
print({"TNF_IL6_NR"})
print({TNF_IL6_NR})
FAP_IL6_R <- cor.test(x=test[which(rownames(test) %in% names(sample.pheno[which(sample.pheno %in% "L_R")])),][,"FAP"],y=test[which(rownames(test) %in% names(sample.pheno[which(sample.pheno %in% "L_R")])),][,"IL6"],method="pearson")
print({"FAP_IL6_R"})
print({FAP_IL6_R})
FAP_IL6_NR <- cor.test(x=test[which(rownames(test) %in% names(sample.pheno[which(sample.pheno %in% "L_NR")])),][,"FAP"],y=test[which(rownames(test) %in% names(sample.pheno[which(sample.pheno %in% "L_NR")])),][,"IL6"],method="pearson")
print({"FAP_IL6_NR"})
print({FAP_IL6_NR})
write.table(TNF_IL6_R[c("estimate","p.value")], paste0(work_dir,"02.Table/",name,"_TNF_IL6_R.txt"), quote=F, sep='\t')
write.table(TNF_IL6_NR[c("estimate","p.value")], paste0(work_dir,"02.Table/",name,"_TNF_IL6_NR.txt"), quote=F, sep='\t')
write.table(FAP_IL6_R[c("estimate","p.value")], paste0(work_dir,"02.Table/",name,"_FAP_IL6_R.txt"), quote=F, sep='\t')
write.table(FAP_IL6_NR[c("estimate","p.value")], paste0(work_dir,"02.Table/",name,"_FAP_IL6_NR.txt"), quote=F, sep='\t')

#FDR and logFC
Ftest_FDR <- topTags(Ftest, n=rownames(Ftest))[,c(7,11)]
Ftest_FDR[which(rownames(Ftest_FDR) %in% Iif_genes),]

#Inflammatory_loops, Il
#Il_genes <- c("IL17","TNF","LIF","LIFR","IFNAR1","IFNAR2","IL6ST","STAT4","JAK1","TYK2","IL6","CSF3","IL33","IL1A","IL1B")
Il_genes <- c("IL6","PDPN","CD90","CXCL8","CXCL1","CXCL5","CXCL2","IL11","CCL2","CCL7","IL33","CD74","CCL19","CCL21")
wIl <- which(y_uniq$genes$gene_name %in% Il_genes)
y_uniq_Il <- y_uniq[wIl,]
y_uniq_Il_heatmap <- y_uniq_Il
rownames(y_uniq_Il_heatmap) <- y_uniq_Il$genes$gene_name
y_uniq_Il.lcpm <- cpm(y_uniq_Il_heatmap, log=T)

#FAP
FAP_genes <- c("FAP","IL6")
wFAP <- which(y_uniq$genes$gene_name %in% FAP_genes)
y_uniq_FAP <- y_uniq[wFAP,]
y_uniq_FAP_heatmap <- y_uniq_FAP
rownames(y_uniq_FAP_heatmap) <- y_uniq_FAP$genes$gene_name
y_uniq_FAP.lcpm <- cpm(y_uniq_FAP_heatmap, log=T)
#FDR and logFC
Ftest_FDR <- topTags(Ftest, n=rownames(Ftest))[,c(7,11)]
Ftest_FDR[which(rownames(Ftest_FDR) %in% FAP_genes),]

#PCA
y_uniq_dt.lcpm.t <- t(y_uniq_dt.lcpm)
y_uniq_dt.pca <- pca(y_uniq_dt.lcpm.t, ncomp=10, center=TRUE, scale=FALSE)
png(filename=paste0(work_dir,"01.Plot/dt_pca.png"), width=20, height=15, units="cm", res=200)
pca.variates <- y_uniq_dt.pca$variates[[1]]
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
write.table(y_uniq_dt.pca$variates, paste0(work_dir,"01.Plot/pca_variates.csv"), quote=F, col.names=T, row.names=T,sep=",")

#Heatmap
y_uniq_dt.lcpm <- t(y_uniq_dt.lcpm)
for(i in 1:length(colnames(y_uniq_dt.lcpm))){
y_uniq_dt.lcpm[,i] <- scale(y_uniq_dt.lcpm[,i], center=TRUE, scale=FALSE)
}
y_uniq_dt.lcpm <- t(y_uniq_dt.lcpm)
cols <- list(pheno=c("L_NR"="blue", "L_R"="orange"))
anno <- HeatmapAnnotation(
    pheno=sample.pheno, col=cols,
    simple_anno_size = unit(1.5, "cm"),height = unit(1.5, "cm"),
    annotation_name_rot = 45,
    annotation_name_gp=gpar(fontsize=20),
    annotation_legend_param=list(
        pheno=list(title_gp=gpar(fontsize=20), labels_gp=gpar(fontsize=20))
    )
)

png(filename=paste0(work_dir,"01.Plot/Heatmap.png"), width=60, height=150, units="cm", res=200)
print({Heatmap(y_uniq_dt.lcpm, top_annotation=anno,
        column_names_gp=grid::gpar(fontsize=12), row_names_gp = grid::gpar(fontsize = 12),
        width=unit(45, "cm"), height=unit(135,"cm"),
        heatmap_legend_param = list(title = "mat", labels_gp = gpar(fontsize = 25), title_gp = gpar(fontsize = 30))
)})
dev.off()

y_uniq_vali.lcpm.scale <- t(y_uniq_vali.lcpm)
#Heatmap
print({paste0("y_uniq_vali.lcpm's dimension is ",dim(y_uniq_vali.lcpm))})
for(i in 1:length(colnames(y_uniq_vali.lcpm.scale))){
y_uniq_vali.lcpm.scale[,i] <- scale(y_uniq_vali.lcpm.scale[,i], center=TRUE, scale=FALSE)
}
y_uniq_vali.lcpm.scale <- t(y_uniq_vali.lcpm.scale)
png(filename=paste0(work_dir,"01.Plot/246_Heatmap.png"), width=50, height=60, units="cm", res=200)
print({Heatmap(y_uniq_vali.lcpm.scale, top_annotation=anno,
        column_names_gp=grid::gpar(fontsize=10), row_names_gp = grid::gpar(fontsize = 10),
        width=unit(42, "cm"), height=unit(50,"cm"),
        heatmap_legend_param = list(title = "mat", labels_gp = gpar(fontsize = 10), title_gp = gpar(fontsize = 10))
)})
dev.off()

#Iif
y_uniq_Iif.lcpm.scale <- t(y_uniq_Iif.lcpm)
#Heatmap
print({paste0("y_uniq_Iif.lcpm's dimension is ",dim(y_uniq_Iif.lcpm))})
for(i in 1:length(colnames(y_uniq_Iif.lcpm.scale))){
y_uniq_Iif.lcpm.scale[,i] <- scale(y_uniq_Iif.lcpm.scale[,i], center=TRUE, scale=FALSE)
}
y_uniq_Iif.lcpm.scale <- t(y_uniq_Iif.lcpm.scale)
png(filename=paste0(work_dir,"01.Plot/Iif_Heatmap.png"), width=50, height=30, units="cm", res=200)
print({Heatmap(y_uniq_Iif.lcpm.scale, top_annotation=anno,
        column_names_gp=grid::gpar(fontsize=10), row_names_gp = grid::gpar(fontsize = 22),
        width=unit(42, "cm"), height=unit(25,"cm"),
        heatmap_legend_param = list(title = "mat", labels_gp = gpar(fontsize = 20), title_gp = gpar(fontsize = 20))
)})
dev.off()

#Il
y_uniq_Il.lcpm.scale <- t(y_uniq_Il.lcpm)
#Heatmap
print({paste0("y_uniq_Il.lcpm's dimension is ",dim(y_uniq_Il.lcpm))})
for(i in 1:length(colnames(y_uniq_Il.lcpm.scale))){
y_uniq_Il.lcpm.scale[,i] <- scale(y_uniq_Il.lcpm.scale[,i], center=TRUE, scale=FALSE)
}
y_uniq_Il.lcpm.scale <- t(y_uniq_Il.lcpm.scale)
png(filename=paste0(work_dir,"01.Plot/Il_Heatmap.png"), width=52, height=30, units="cm", res=200)
print({Heatmap(y_uniq_Il.lcpm.scale, top_annotation=anno,
        column_names_gp=grid::gpar(fontsize=10), row_names_gp = grid::gpar(fontsize = 22),
        width=unit(42, "cm"), height=unit(25,"cm"),
        heatmap_legend_param = list(title = "mat", labels_gp = gpar(fontsize = 10), title_gp = gpar(fontsize = 10))
)})
dev.off()




}

