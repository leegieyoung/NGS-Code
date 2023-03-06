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
#y_uniq <- estimateDisp(y_uniq,sample.design)
y_uniq <- estimateDisp(y_uniq) #Using classic mode
print({
	paste0("Do Norm is :" ,dim(y_uniq))
	})
y_uniq_cpm <- cpm(y_uniq, log=F)
write.table(y_uniq_cpm, paste0(work_dir,"02.Table/",name,"_cpm.txt"), quote=F, sep='\t')

#####Analysis####
print("Run glmQy_uniqLF")
#715 genes table
test <- read.table("/data/keeyoung/gy_RNA/01.Samplelist/715genes.txt", sep='\t')
colnames(test) <- test[1,]
test <- test[-1,]
Total_genes <- test$Gene
iCMS2_Up <- test$Gene[which(test$Type %in% "iCMS2_Up")]
iCMS3_Down <- test$Gene[which(test$Type %in% "iCMS3_Down")]
iCMS2_Down <- test$Gene[which(test$Type %in% "iCMS2_Down")]
iCMS3_Up <- test$Gene[which(test$Type %in% "iCMS3_Up")]

###heatmap############
#cols <- list(pheno=c("iCMS2"="blue", "iCMS3"="pink", "Intermediate"="gray"))
#anno <- HeatmapAnnotation(
#    pheno=sample.pheno, col=cols,
#    simple_anno_size = unit(1.5, "cm"),height = unit(1.5, "cm"),
#    annotation_name_rot = 45,
#    annotation_name_gp=gpar(fontsize=20),
#    annotation_legend_param=list(
#        pheno=list(title_gp=gpar(fontsize=20), labels_gp=gpar(fontsize=20))
#    )
#)

#cols <- list(iCMS=c("iCMS2"="blue", "iCMS3"="pink"), MSI=c("Stability"="green", "Unknwon"="gray","Instability"="yellow"), Type=c("Tumor"="Red", "Adenoma"="orange"))
#anno <- HeatmapAnnotation(iCMS=(iCMS=sample.pheno), CMS2=anno_barplot(sample.CMS2), CMS3=anno_barplot(sample.CMS3), CMS1=anno_barplot(sample.CMS1), MSI=(MSI=sample.MSI), Type=(Type=sample.Type), col=cols,
#	simple_anno_size = unit(1.5, "cm"),height = unit(15, "cm"), annotation_name_gp=gpar(fontsize=35),
#	annotation_legend_param=list(iCMS=list(title_gp=gpar(fontsize=35), labels_gp=gpar(fontsize=30)), MSI=list(title_gp=gpar(fontsize=35), labels_gp=gpar(fontsize=30)), Type=list(title_gp=gpar(fontsize=35), labels_gp=gpar(fontsize=30)))
#)
cols <- list(iCMS=c("iCMS2"="blue", "iCMS3"="pink"), CMS_frac=c("blue","pink","yellow","green"),
	MS=c("Stability"="black", "Unknwon"="gray","Instability"="yellow"), Type=c("Tumor"="Red", "Adenoma"="orange")
)
anno <- HeatmapAnnotation(iCMS=(iCMS=sample.pheno), CMS_frac=anno_barplot(cbind(sample.CMS2,sample.CMS3,sample.CMS1,sample.CMS4), gp=gpar(fill=c("blue","pink","yellow","green"))),
	MS=(MS=sample.MSI), Type=(Type=sample.Type), col=cols,
    simple_anno_size = unit(2.5, "cm"),height = unit(20, "cm"), annotation_name_gp=gpar(fontsize=60, fontface="bold"),
    annotation_legend_param=list(iCMS=list(title_gp=gpar(fontsize=55, fontface="bold"), labels_gp=gpar(fontsize=50)), MS=list(title_gp=gpar(fontsize=55, fontface="bold"), labels_gp=gpar(fontsize=50)), Type=list(title_gp=gpar(fontsize=55, fontface="bold"), labels_gp=gpar(fontsize=50)))
)

#The 715 order genes
w_vali <- which(y_uniq$genes$gene_name %in% Total_genes)
y_uniq_vali <- y_uniq[w_vali,]
y_uniq_vali_heatmap <- y_uniq_vali
rownames(y_uniq_vali_heatmap) <- y_uniq_vali$genes$gene_name
y_uniq_vali.lcpm <- cpm(y_uniq_vali_heatmap, log=T)
colnames(y_uniq_vali.lcpm) <- gsub(".TR","",colnames(y_uniq_vali.lcpm))
y_uniq_vali.lcpm.scale <- t(y_uniq_vali.lcpm)
print({paste0("y_uniq_vali.lcpm's dimension is ",dim(y_uniq_vali.lcpm))})
for(i in 1:length(colnames(y_uniq_vali.lcpm.scale))){
y_uniq_vali.lcpm.scale[,i] <- scale(y_uniq_vali.lcpm.scale[,i], center=TRUE, scale=FALSE)
}
y_uniq_vali.lcpm.scale <- t(y_uniq_vali.lcpm.scale)
#png(filename=paste0(work_dir,"01.Plot/715_Heatmap.png"), width=50, height=60, units="cm", res=200)
Total_genes <- Total_genes[which(Total_genes %in% rownames(y_uniq_vali.lcpm.scale))]
pdf(paste0(work_dir,"01.Plot/Total_genes_Heatmap.pdf"), width=25, height=30)
print({Heatmap(y_uniq_vali.lcpm.scale, top_annotation=anno,
        column_names_gp=grid::gpar(fontsize=25), row_names_gp = grid::gpar(fontsize = 0.1),
		row_order = Total_genes,
        width=unit(42, "cm"), height=unit(25,"cm"),
        heatmap_legend_param = list(title = "Mat", labels_gp = gpar(fontsize = 25), title_gp = gpar(fontsize = 35), border="red",grid_width = unit(1, "cm")),
		column_title = "715_ordered genes(used 674 genes)", column_title_gp = gpar(fontsize = 40, fontface = "bold")
)})
dev.off()

pathway.gene <- scan("/data/keeyoung/gy_RNA/Code/Oncoprint/pathway.gene", what=character(0))

crc <- read.table("/data/keeyoung/gy_RNA/CRC/05.variant_calling/26df2.txt", head=TRUE,stringsAsFactors=FALSE,sep='\t')
rownames(crc) <- crc[,1]
crc <- crc[,-1]
all.equal(colnames(crc[,match(colnames(y_uniq_vali.lcpm), colnames(crc))]) , colnames(y_uniq_vali.lcpm))

crc <- crc[,match(colnames(y_uniq_vali.lcpm), colnames(crc))]
# just for demonstration
oncol = c("nonsense" = "red", "nonsynonymous" = "blue", "frameshift" = "#008000")
alter_fun = list(
    background = alter_graphic("rect", width = 0.95, height = 0.75, fill = "#CCCCCC"),   
    nonsense = alter_graphic("rect", width = 0.8,height=0.75, fill = oncol["nonsense"]),
    nonsynonymous = alter_graphic("rect", width = 0.33,height=0.75, fill = oncol["nonsynonymous"]),
    frameshift = alter_graphic("rect", width = 0.8, height = 0.33, fill = oncol["frameshift"])
)
heatmap_legend_param = list(title = "Alternations", at = c("nonsense", "nonsynonymous", "frameshift"),  
        labels = c("Nonsense", "Nonsynonymous", "Frameshift"), 
		labels_gp = gpar(fontsize = 45), title_gp = gpar(fontsize = 55, fontface="bold"), grid_width = unit(3, "cm"))

ht_list <- 
	oncoPrint(crc,
    	alter_fun = alter_fun, col = oncol, column_order=gsub(".TR","",names(sample.pheno)), row_names_gp = grid::gpar(fontsize = 45, "cm"),
	    column_title = "715_ordered genes(used 674 genes)", column_title_gp = gpar(fontsize = 75, fontface = "bold"), row_order=pathway.gene,
		top_annotation = HeatmapAnnotation(column_barplot = anno_oncoprint_barplot(border = TRUE, height = unit(8, "cm"),show_fraction = FALSE)),
		right_annotation = rowAnnotation(row_barplot = anno_oncoprint_barplot(border = TRUE, width = unit(8, "cm"), show_fraction = FALSE, ,axis_param = list(side = "bottom", labels_rot = 90))),
		heatmap_legend_param = heatmap_legend_param, alter_fun_is_vectorized = FALSE, height=unit(35,"cm"),
		pct_gp= grid::gpar(fontsize = 45, "cm")) %v%
	Heatmap(y_uniq_vali.lcpm.scale, top_annotation=anno,
        column_names_gp=grid::gpar(fontsize=55), row_names_gp = grid::gpar(fontsize = 0.1),
        row_order = Total_genes,
        width=unit(25, "cm"), height=unit(30,"cm"),
        heatmap_legend_param = list(title = "Mat", labels_gp = gpar(fontsize = 45), title_gp = gpar(fontsize = 55, fontface="bold"), border="red",grid_width = unit(3, "cm")),
		)
pdf(paste0(work_dir,"01.Plot/Total_genes_Heatmap_Oncoprint.pdf"), width=60, height=60)
print({draw(ht_list, merge_legend=TRUE, heatmap_legend_side = "left", annotation_legend_side = "bottom")})
#print({ht_list})
dev.off()
pdf(paste0(work_dir,"01.Plot/Total_genes_Heatmap.pdf"), width=25, height=30)

########################################################
#############Article##########################################
oncol = c("nonsense" = "red", "nonsynonymous" = "blue", "frameshift" = "#008000")
alter_fun = list(
    background = alter_graphic("rect", width = 0.65, height = 0.75, fill = "#CCCCCC"),
    nonsense = alter_graphic("rect", width = 0.65, ,height=0.75, fill = oncol["nonsense"]),
    nonsynonymous = alter_graphic("rect", width = 0.28,height=0.75, fill = oncol["nonsynonymous"]),
    frameshift = alter_graphic("rect", width = 0.65, height = 0.33, fill = oncol["frameshift"])
)
ht_list <-
    oncoPrint(crc,
        alter_fun = alter_fun, col = oncol, column_order=gsub(".TR","",names(sample.pheno)), row_names_gp = grid::gpar(fontsize = 45, "cm"),
        column_title = "715_ordered genes(used 674 genes)", column_title_gp = gpar(fontsize = 75, fontface = "bold"), row_order=pathway.gene,
        top_annotation = HeatmapAnnotation(column_barplot = anno_oncoprint_barplot(border = TRUE, height = unit(8, "cm"),show_fraction = FALSE)),
        right_annotation = rowAnnotation(row_barplot = anno_oncoprint_barplot(border = TRUE, width = unit(8, "cm"), show_fraction = FALSE, ,axis_param = list(side = "bottom", labels_rot = 90))),
        heatmap_legend_param = heatmap_legend_param, alter_fun_is_vectorized = FALSE, height=unit(35,"cm"),
        pct_gp= grid::gpar(fontsize = 45, "cm")) %v%
    Heatmap(y_uniq_vali.lcpm.scale, top_annotation=anno,
		show_column_names=FALSE, show_row_names=FALSE,
        row_order = Total_genes,
        width=unit(15, "cm"), height=unit(30,"cm"),
        heatmap_legend_param = list(title = "Mat", labels_gp = gpar(fontsize = 45), title_gp = gpar(fontsize = 55, fontface="bold"), border="red",grid_width = unit(3, "cm")),
        )

pdf(paste0(work_dir,"01.Plot/Total_genes_Heatmap_Oncoprint_article.pdf"), width=30, height=60)
print({draw(ht_list, merge_legend=TRUE, heatmap_legend_side = "left", annotation_legend_side = "bottom")})
#print({ht_list})
dev.off()

###############################################################
#iCMS2
iCMS2.pheno <- sample.pheno[which(sample.pheno %in% "iCMS2")]
names(iCMS2.pheno) <- gsub(".TR","",names(iCMS2.pheno))
crc.iCMS2 <- crc[,which(colnames(crc) %in% names(iCMS2.pheno))]
y_uniq_vali.lcpm.scale.iCMS2 <- y_uniq_vali.lcpm.scale[,which(colnames(y_uniq_vali.lcpm.scale) %in% names(iCMS2.pheno))]

ht_list_iCMS2 <-
    oncoPrint(crc.iCMS2,
        alter_fun = alter_fun, col = oncol, 
		#column_order=gsub(".TR","",names(iCMS2.pheno)), 
		row_names_gp = grid::gpar(fontsize = 25, "cm"),
        column_title = "The only iCMS2", column_title_gp = gpar(fontsize = 60, fontface = "bold"), row_order=pathway.gene,
        heatmap_legend_param = heatmap_legend_param, alter_fun_is_vectorized = FALSE, height=unit(35,"cm")) %v%
    Heatmap(y_uniq_vali.lcpm.scale.iCMS2,
        column_names_gp=grid::gpar(fontsize=25), row_names_gp = grid::gpar(fontsize = 0.1),
        row_order = Total_genes,
        width=unit(12, "cm"), height=unit(25,"cm"),
        heatmap_legend_param = list(title = "Mat", labels_gp = gpar(fontsize = 25), title_gp = gpar(fontsize = 35), border="red",grid_width = unit(1, "cm")),
        )
pdf(paste0(work_dir,"01.Plot/Total_genes_Heatmap_Oncoprint_iCMS2.pdf"), width=40, height=60)
print({ht_list_iCMS2})
dev.off()

#iCMS3
iCMS3.pheno <- sample.pheno[which(sample.pheno %in% "iCMS3")]
names(iCMS3.pheno) <- gsub(".TR","",names(iCMS3.pheno))
crc.iCMS3 <- crc[,which(colnames(crc) %in% names(iCMS3.pheno))]
y_uniq_vali.lcpm.scale.iCMS3 <- y_uniq_vali.lcpm.scale[,which(colnames(y_uniq_vali.lcpm.scale) %in% names(iCMS3.pheno))]

ht_list_iCMS3 <-
    oncoPrint(crc.iCMS3,
        alter_fun = alter_fun, col = oncol, 
        #column_order=gsub(".TR","",names(iCMS3.pheno)), 
        row_names_gp = grid::gpar(fontsize = 25, "cm"),
        column_title = "The only iCMS3", column_title_gp = gpar(fontsize = 60, fontface = "bold"), row_order=pathway.gene,
        heatmap_legend_param = heatmap_legend_param, alter_fun_is_vectorized = FALSE, height=unit(35,"cm")) %v%
    Heatmap(y_uniq_vali.lcpm.scale.iCMS3,
        column_names_gp=grid::gpar(fontsize=25), row_names_gp = grid::gpar(fontsize = 0.1),
        row_order = Total_genes,
        width=unit(12, "cm"), height=unit(25,"cm"),
        heatmap_legend_param = list(title = "Mat", labels_gp = gpar(fontsize = 25), title_gp = gpar(fontsize = 35), border="red",grid_width = unit(1, "cm")),
        )  
pdf(paste0(work_dir,"01.Plot/Total_genes_Heatmap_Oncoprint_iCMS3.pdf"), width=40, height=60)
print({ht_list_iCMS3})
dev.off()

#box plot
myo <- read.table("/data/keeyoung/gy_RNA/CRC/Cibersortx_myofibroblasts/stromal_Myofibroblasts.txt", sep='\t', head=T)
rownames(myo) <- myo[,1]
myo <- myo[,-1]

rmna.myo <- na.omit(myo)
p <- ggplot(myo, aes(x=RF.nearestCMS, y=Myofibroblasts)) + geom_boxplot()
q <- ggplot(rmna.myo, aes(x=RF.predictedCMS, y=Myofibroblasts)) + geom_boxplot()
pdf(paste0(work_dir,"01.Plot/boxplot_nearest.pdf"), width=25, height=30)
print({p})
dev.off()

pdf(paste0(work_dir,"01.Plot/boxplot_predicted.pdf"), width=25, height=30)
print({q})
dev.off()


######################################################################################old#################################
#The iCMS2_Up of 715 Genes 
w_vali <- which(y_uniq$genes$gene_name %in% iCMS2_Up)
y_uniq_vali <- y_uniq[w_vali,]
y_uniq_vali_heatmap <- y_uniq_vali
rownames(y_uniq_vali_heatmap) <- y_uniq_vali$genes$gene_name
y_uniq_vali.lcpm <- cpm(y_uniq_vali_heatmap, log=T)

y_uniq_vali.lcpm.scale <- t(y_uniq_vali.lcpm)
print({paste0("y_uniq_vali.lcpm's dimension is ",dim(y_uniq_vali.lcpm))})
for(i in 1:length(colnames(y_uniq_vali.lcpm.scale))){
y_uniq_vali.lcpm.scale[,i] <- scale(y_uniq_vali.lcpm.scale[,i], center=TRUE, scale=FALSE)
}
y_uniq_vali.lcpm.scale <- t(y_uniq_vali.lcpm.scale)
#png(filename=paste0(work_dir,"01.Plot/715_Heatmap.png"), width=50, height=60, units="cm", res=200)
pdf(paste0(work_dir,"01.Plot/iCMS2_Up_Heatmap.pdf"), width=25, height=30)
print({Heatmap(y_uniq_vali.lcpm.scale, top_annotation=anno,
        column_names_gp=grid::gpar(fontsize=10), row_names_gp = grid::gpar(fontsize = 10),
        width=unit(42, "cm"), height=unit(50,"cm"),
        heatmap_legend_param = list(title = "mat", labels_gp = gpar(fontsize = 10), title_gp = gpar(fontsize = 10)),
		column_title = "iCMS2_Up", column_title_gp = gpar(fontsize = 40, fontface = "bold")
)})
dev.off()

#The iCMS3_Down of 715 Genes 
w_vali <- which(y_uniq$genes$gene_name %in% iCMS3_Down)
y_uniq_vali <- y_uniq[w_vali,]
y_uniq_vali_heatmap <- y_uniq_vali
rownames(y_uniq_vali_heatmap) <- y_uniq_vali$genes$gene_name
y_uniq_vali.lcpm <- cpm(y_uniq_vali_heatmap, log=T)

y_uniq_vali.lcpm.scale <- t(y_uniq_vali.lcpm)
print({paste0("y_uniq_vali.lcpm's dimension is ",dim(y_uniq_vali.lcpm))})
for(i in 1:length(colnames(y_uniq_vali.lcpm.scale))){
y_uniq_vali.lcpm.scale[,i] <- scale(y_uniq_vali.lcpm.scale[,i], center=TRUE, scale=FALSE)
}
y_uniq_vali.lcpm.scale <- t(y_uniq_vali.lcpm.scale)
#png(filename=paste0(work_dir,"01.Plot/715_Heatmap.png"), width=50, height=60, units="cm", res=200)
pdf(paste0(work_dir,"01.Plot/iCMS3_Down_Heatmap.pdf"), width=25, height=30)
print({Heatmap(y_uniq_vali.lcpm.scale, top_annotation=anno,
        column_names_gp=grid::gpar(fontsize=10), row_names_gp = grid::gpar(fontsize = 10),
        width=unit(42, "cm"), height=unit(50,"cm"),
        heatmap_legend_param = list(title = "mat", labels_gp = gpar(fontsize = 10), title_gp = gpar(fontsize = 10)),
        column_title = "iCMS3_Down", column_title_gp = gpar(fontsize = 40, fontface = "bold")
)})
dev.off()

#The iCMS2_Down of 715 Genes 
w_vali <- which(y_uniq$genes$gene_name %in% iCMS2_Down)
y_uniq_vali <- y_uniq[w_vali,]
y_uniq_vali_heatmap <- y_uniq_vali
rownames(y_uniq_vali_heatmap) <- y_uniq_vali$genes$gene_name
y_uniq_vali.lcpm <- cpm(y_uniq_vali_heatmap, log=T)

y_uniq_vali.lcpm.scale <- t(y_uniq_vali.lcpm)
print({paste0("y_uniq_vali.lcpm's dimension is ",dim(y_uniq_vali.lcpm))})
for(i in 1:length(colnames(y_uniq_vali.lcpm.scale))){
y_uniq_vali.lcpm.scale[,i] <- scale(y_uniq_vali.lcpm.scale[,i], center=TRUE, scale=FALSE)
}
y_uniq_vali.lcpm.scale <- t(y_uniq_vali.lcpm.scale)
#png(filename=paste0(work_dir,"01.Plot/715_Heatmap.png"), width=50, height=60, units="cm", res=200)
pdf(paste0(work_dir,"01.Plot/iCMS2_Down_Heatmap.pdf"), width=25, height=30)
print({Heatmap(y_uniq_vali.lcpm.scale, top_annotation=anno,
        column_names_gp=grid::gpar(fontsize=10), row_names_gp = grid::gpar(fontsize = 10),
        width=unit(42, "cm"), height=unit(50,"cm"),
        heatmap_legend_param = list(title = "mat", labels_gp = gpar(fontsize = 10), title_gp = gpar(fontsize = 10)),
        column_title = "iCMS2_Down", column_title_gp = gpar(fontsize = 40, fontface = "bold")
)})
dev.off()

#The iCMS3_Up of 715 Genes 
w_vali <- which(y_uniq$genes$gene_name %in% iCMS3_Up)
y_uniq_vali <- y_uniq[w_vali,]
y_uniq_vali_heatmap <- y_uniq_vali
rownames(y_uniq_vali_heatmap) <- y_uniq_vali$genes$gene_name
y_uniq_vali.lcpm <- cpm(y_uniq_vali_heatmap, log=T)

y_uniq_vali.lcpm.scale <- t(y_uniq_vali.lcpm)
print({paste0("y_uniq_vali.lcpm's dimension is ",dim(y_uniq_vali.lcpm))})
for(i in 1:length(colnames(y_uniq_vali.lcpm.scale))){
y_uniq_vali.lcpm.scale[,i] <- scale(y_uniq_vali.lcpm.scale[,i], center=TRUE, scale=FALSE)
}
y_uniq_vali.lcpm.scale <- t(y_uniq_vali.lcpm.scale)
#png(filename=paste0(work_dir,"01.Plot/715_Heatmap.png"), width=50, height=60, units="cm", res=200)
pdf(paste0(work_dir,"01.Plot/iCMS3_Up_Heatmap.pdf"), width=25, height=30)
print({Heatmap(y_uniq_vali.lcpm.scale, top_annotation=anno,
        column_names_gp=grid::gpar(fontsize=10), row_names_gp = grid::gpar(fontsize = 10),
        width=unit(42, "cm"), height=unit(50,"cm"),
        heatmap_legend_param = list(title = "mat", labels_gp = gpar(fontsize = 10), title_gp = gpar(fontsize = 10)),
        column_title = "iCMS3_Up", column_title_gp = gpar(fontsize = 40, fontface = "bold")
)})
dev.off()

sample.order <- scan("/data/keeyoung/gy_RNA/01.Samplelist/CRCsamples.order", what=character(0))
#The iCMS2_Up of 715 Genes 
w_vali <- which(y_uniq$genes$gene_name %in% iCMS2_Up)
y_uniq_vali <- y_uniq[w_vali,]
y_uniq_vali_heatmap <- y_uniq_vali
rownames(y_uniq_vali_heatmap) <- y_uniq_vali$genes$gene_name
y_uniq_vali.lcpm <- cpm(y_uniq_vali_heatmap, log=T)

y_uniq_vali.lcpm.scale <- t(y_uniq_vali.lcpm)
print({paste0("y_uniq_vali.lcpm's dimension is ",dim(y_uniq_vali.lcpm))})
for(i in 1:length(colnames(y_uniq_vali.lcpm.scale))){
y_uniq_vali.lcpm.scale[,i] <- scale(y_uniq_vali.lcpm.scale[,i], center=TRUE, scale=FALSE)
}
y_uniq_vali.lcpm.scale <- t(y_uniq_vali.lcpm.scale)
#png(filename=paste0(work_dir,"01.Plot/715_Heatmap.png"), width=50, height=60, units="cm", res=200)
pdf(paste0(work_dir,"01.Plot/iCMS2_Up_Heatmap_order.pdf"), width=25, height=30)
print({Heatmap(y_uniq_vali.lcpm.scale, top_annotation=anno,
		column_order = sample.order,
        column_names_gp=grid::gpar(fontsize=10), row_names_gp = grid::gpar(fontsize = 10),
        width=unit(42, "cm"), height=unit(50,"cm"),
        heatmap_legend_param = list(title = "mat", labels_gp = gpar(fontsize = 10), title_gp = gpar(fontsize = 10)),
        column_title = "Reorder Pheno, Genes are iCMS2_Up", column_title_gp = gpar(fontsize = 40, fontface = "bold")
)})
dev.off()

#The iCMS3_Down of 715 Genes 
w_vali <- which(y_uniq$genes$gene_name %in% iCMS3_Down)
y_uniq_vali <- y_uniq[w_vali,]
y_uniq_vali_heatmap <- y_uniq_vali
rownames(y_uniq_vali_heatmap) <- y_uniq_vali$genes$gene_name
y_uniq_vali.lcpm <- cpm(y_uniq_vali_heatmap, log=T)

y_uniq_vali.lcpm.scale <- t(y_uniq_vali.lcpm)
print({paste0("y_uniq_vali.lcpm's dimension is ",dim(y_uniq_vali.lcpm))})
for(i in 1:length(colnames(y_uniq_vali.lcpm.scale))){
y_uniq_vali.lcpm.scale[,i] <- scale(y_uniq_vali.lcpm.scale[,i], center=TRUE, scale=FALSE)
}
y_uniq_vali.lcpm.scale <- t(y_uniq_vali.lcpm.scale)
#png(filename=paste0(work_dir,"01.Plot/715_Heatmap.png"), width=50, height=60, units="cm", res=200)
pdf(paste0(work_dir,"01.Plot/iCMS3_Down_Heatmap_order.pdf"), width=25, height=30)
print({Heatmap(y_uniq_vali.lcpm.scale, top_annotation=anno,
        column_order = sample.order,
        column_names_gp=grid::gpar(fontsize=10), row_names_gp = grid::gpar(fontsize = 10),
        width=unit(42, "cm"), height=unit(50,"cm"),
        heatmap_legend_param = list(title = "mat", labels_gp = gpar(fontsize = 10), title_gp = gpar(fontsize = 10)),
        column_title = "Reorder Pheno, Genes are iCMS3_Down", column_title_gp = gpar(fontsize = 40, fontface = "bold")
)})  
dev.off()

#The iCMS2_Down of 715 Genes 
w_vali <- which(y_uniq$genes$gene_name %in% iCMS2_Down)
y_uniq_vali <- y_uniq[w_vali,]
y_uniq_vali_heatmap <- y_uniq_vali
rownames(y_uniq_vali_heatmap) <- y_uniq_vali$genes$gene_name
y_uniq_vali.lcpm <- cpm(y_uniq_vali_heatmap, log=T)

y_uniq_vali.lcpm.scale <- t(y_uniq_vali.lcpm)
print({paste0("y_uniq_vali.lcpm's dimension is ",dim(y_uniq_vali.lcpm))})
for(i in 1:length(colnames(y_uniq_vali.lcpm.scale))){
y_uniq_vali.lcpm.scale[,i] <- scale(y_uniq_vali.lcpm.scale[,i], center=TRUE, scale=FALSE)
}
y_uniq_vali.lcpm.scale <- t(y_uniq_vali.lcpm.scale)
#png(filename=paste0(work_dir,"01.Plot/715_Heatmap.png"), width=50, height=60, units="cm", res=200)
pdf(paste0(work_dir,"01.Plot/iCMS2_Down_Heatmap_order.pdf"), width=25, height=30)
print({Heatmap(y_uniq_vali.lcpm.scale, top_annotation=anno,
        column_order = sample.order,
        column_names_gp=grid::gpar(fontsize=10), row_names_gp = grid::gpar(fontsize = 10),
        width=unit(42, "cm"), height=unit(50,"cm"),
        heatmap_legend_param = list(title = "mat", labels_gp = gpar(fontsize = 10), title_gp = gpar(fontsize = 10)),
        column_title = "Reorder Pheno, Genes are iCMS2_Down", column_title_gp = gpar(fontsize = 40, fontface = "bold")
)})  
dev.off()

#The iCMS3_Up of 715 Genes 
w_vali <- which(y_uniq$genes$gene_name %in% iCMS3_Up)
y_uniq_vali <- y_uniq[w_vali,]
y_uniq_vali_heatmap <- y_uniq_vali
rownames(y_uniq_vali_heatmap) <- y_uniq_vali$genes$gene_name
y_uniq_vali.lcpm <- cpm(y_uniq_vali_heatmap, log=T)

y_uniq_vali.lcpm.scale <- t(y_uniq_vali.lcpm)
print({paste0("y_uniq_vali.lcpm's dimension is ",dim(y_uniq_vali.lcpm))})
for(i in 1:length(colnames(y_uniq_vali.lcpm.scale))){
y_uniq_vali.lcpm.scale[,i] <- scale(y_uniq_vali.lcpm.scale[,i], center=TRUE, scale=FALSE)
}
y_uniq_vali.lcpm.scale <- t(y_uniq_vali.lcpm.scale)
#png(filename=paste0(work_dir,"01.Plot/715_Heatmap.png"), width=50, height=60, units="cm", res=200)
pdf(paste0(work_dir,"01.Plot/iCMS3_Up_Heatmap_order.pdf"), width=25, height=30)
print({Heatmap(y_uniq_vali.lcpm.scale, top_annotation=anno,
        column_order = sample.order,
        column_names_gp=grid::gpar(fontsize=10), row_names_gp = grid::gpar(fontsize = 10),
        width=unit(42, "cm"), height=unit(50,"cm"),
        heatmap_legend_param = list(title = "mat", labels_gp = gpar(fontsize = 10), title_gp = gpar(fontsize = 10)),
        column_title = "Reorder Pheno, Genes are iCMS3_Up", column_title_gp = gpar(fontsize = 40, fontface = "bold")
)})  
dev.off()


}

