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
library(gridExtra)
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

cols <- list(iCMS=c("iCMS2"="#0040FF", "iCMS3"="#F5A9D0"), CMS.frac=c("yellow","#0040FF","#F5A9D0","green"),
    MSI.status=c("Stability"="#A9A9F5", "Unknwon"="gray","Instability"="#F2F5A9"), Tissue.Type=c("Tumor"="Red", "Adenoma"="#FACC2E"))

anno <- HeatmapAnnotation(iCMS=(iCMS=sample.pheno), CMS.frac=anno_barplot(cbind(sample.CMS2,sample.CMS3,sample.CMS1,sample.CMS4), gp=gpar(fill=c("blue","pink","yellow","green")), axis_param=list(gp=gpar(fontsize = 20))),
    MSI.status=(MSI=sample.MSI), Tissue.Type=(Type=sample.Type),  col=cols,
    simple_anno_size = unit(1.0, "cm"),height = unit(15, "cm"), annotation_name_gp=gpar(fontsize=25, "cm"),
    annotation_legend_param=list(
		iCMS=list(title_gp=gpar(fontsize=25), labels_gp=gpar(fontsize=22), grid_width=unit(1, "cm")), 
		MSI.status=list(title_gp=gpar(fontsize=25), labels_gp=gpar(fontsize=22), grid_width=unit(1, "cm")), 
		Tissue.Type=list(title_gp=gpar(fontsize=25), labels_gp=gpar(fontsize=22), grid_width=unit(1, "cm"))),
		gap=unit(4,"points")
		#CMS.frac=list(title_gp=gpar(fontsize=25), labels_gp=gpar(fontsize=25), grid_width=unit(1, "cm")))
		#CMS.frac=list(c("CMS2","CMS3","CMS1","CMS4"), title="CMS.frac", legend_gp=gpar(fill=c("#0040FF","#F5A9D0","yellow","green"))))
)
#pdf(paste0(work_dir,"01.Plot/Total_genes_Heatmap_230215.pdf"), width=25, height=15)
pdf(paste0(work_dir,"01.Plot/edit2)Total_genes_Heatmap_230215.pdf"), width=25, height=15)
figure2 <-Heatmap(y_uniq_vali.lcpm.scale, top_annotation=anno, 
		#show_fraction=FALSE,
		#row_names_gp= grid::gpar(fontsize = 25, "cm"),
        column_names_gp=grid::gpar(fontsize=25), row_names_gp = grid::gpar(fontsize = 0.1),
		row_order = Total_genes,
        width=unit(42, "cm"), height=unit(15,"cm"),
		heatmap_legend_param = list(title = "expr", at=c(-4,0,4), labels=c("low","zero","high"), legend_height=unit(4,"cm"), labels_gp = gpar(fontsize = 22), title_gp = gpar(fontsize = 25),grid_width = unit(1.0, "cm"), title_position = "lefttop-rot")
		#column_title = "715_ordered genes(used 674 genes)", column_title_gp = gpar(fontsize = 45, fontface = "bold")
)
draw(figure2, heatmap_legend_list=Legend(title="CMS.frac", labels=c("CMS1","CMS2","CMS3","CMS4"), legend_gp=gpar(fill=c("yellow","#0040FF","#F5A9D0","green")), labels_gp = gpar(fontsize = 22), title_gp = gpar(fontsize = 25),grid_width = unit(1.0, "cm")))
dev.off()

pathway.gene <- scan("/data/keeyoung/gy_RNA/Code/Oncoprint/pathway.gene", what=character(0))
crc <- read.table("/data/keeyoung/gy_RNA/CRC/05.variant_calling/26df2.txt", head=TRUE,stringsAsFactors=FALSE,sep='\t')
rownames(crc) <- crc[,1]
crc <- crc[,-1]
all.equal(colnames(crc[,match(colnames(y_uniq_vali.lcpm), colnames(crc))]) , colnames(y_uniq_vali.lcpm))
crc <- crc[,match(colnames(y_uniq_vali.lcpm), colnames(crc))]

crc <- crc[match(pathway.gene, rownames(crc)),]
all.equal(pathway.gene, rownames(crc))

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
dev.off()

#230207 Oncoprint of CRC article figure
#crc <- read.table("/data/keeyoung/gy_RNA/CRC/05.variant_calling/pathway_nonames_25df2.txt", head=TRUE,stringsAsFactors=FALSE,sep='\t')
crc <- read.table("/data/keeyoung/gy_RNA/CRC/05.variant_calling/pathway_nonames_26df2_tp53.txt", head=TRUE,stringsAsFactors=FALSE,sep='\t') #230303 tp53
rownames(crc) <- crc[,1]
crc <- crc[,-1]
all.equal(colnames(crc[,match(colnames(y_uniq_vali.lcpm), colnames(crc))]) , colnames(y_uniq_vali.lcpm))
crc <- crc[,match(colnames(y_uniq_vali.lcpm), colnames(crc))]

MSI.MSI <- sort(sample.MSI)
MSI.Type <- sample.Type[match(names(MSI.MSI), names(sample.Type))]
all.equal(names(MSI.MSI), names(MSI.Type))
MSI.iCMS <- sample.pheno[match(names(MSI.MSI), names(sample.pheno))]
all.equal(names(MSI.MSI), names(MSI.iCMS))
MSI.column <- gsub(".TR","",names(MSI.MSI))

MSI.crc <- crc[,match(MSI.column, colnames(crc))]

WNT.gene <- scan("/data/keeyoung/gy_RNA/Code/Oncoprint/WNT.gene", what=character(0))
MAPK.gene <- scan("/data/keeyoung/gy_RNA/Code/Oncoprint/MAPK.gene", what=character(0))
TGFB.gene <- scan("/data/keeyoung/gy_RNA/Code/Oncoprint/TGFB.gene", what=character(0))
TP53.gene <- "TP53"

#rep_pathway.gene = scan("/data/keeyoung/gy_RNA/Code/Oncoprint/rep_pathway.gene", what=character(0)) 
oncol = c("nonsense" = "#DF0174", "nonsynonymous" = "#2E9AFE", "frameshift" = "#008000")
alter_fun = list(
    background = alter_graphic("rect", width = 0.85, height = 0.85, fill = "#CCCCCC"),
    nonsense = alter_graphic("rect", width = 0.8, ,height=0.8, fill = oncol["nonsense"]),
    nonsynonymous = alter_graphic("rect", width = 0.28, height=0.8, fill = oncol["nonsynonymous"]),
    frameshift = alter_graphic("rect", width = 0.8, height = 0.28, fill = oncol["frameshift"])
)

cols <- list(iCMS=c("iCMS2"="#0040FF", "iCMS3"="#F5A9D0"), CMS_frac=c("#0040FF","#F5A9D0","yellow","green"),
    MSI.status=c("Stability"="#A9A9F5", "Unknwon"="gray","Instability"="#F2F5A9"), Tissue.Type=c("Tumor"="Red", "Adenoma"="#FACC2E")
)
heatmap_legend_param = list(title = "Alternations", at = c("nonsense", "nonsynonymous", "frameshift"),  
        labels = c("Nonsense", "Nonsynonymous", "Frameshift"), 
        labels_gp = gpar(fontsize = 25), title_gp = gpar(fontsize = 25), grid_width = unit(3, "cm"))

ht_list <- oncoPrint(MSI.crc,
        alter_fun = alter_fun, 
		alter_fun_is_vectorized = FALSE,
		col = oncol,
        column_order = gsub(".TR","",names(MSI.MSI)),
        row_title_gp = grid::gpar(fontsize = 35, "cm"), #Left Gene font size
		row_names_gp = grid::gpar(fontsize = 25, "cm"), #Right Gene font size
		remove_empty_rows = FALSE,
        row_order = rownames(MSI.crc),
#		row_split=c(rep("WNT", length(WNT.gene)),rep("MAPK", length(MAPK.gene)), rep("TGFB", length(TGFB.gene))),
		row_split=c(rep("WNT", length(WNT.gene)),rep("MAPK", length(MAPK.gene)), rep("TGFB", length(TGFB.gene)), rep("TP53", length(TP53.gene))),
		row_gap = unit(10,"mm"),
#		remove_empty_rows = TRUE,
        column_title = "Oncoprint of CRC", column_title_gp = gpar(fontsize = 35, fontface = "bold"), 
		top_annotation = HeatmapAnnotation(column_barplot = anno_oncoprint_barplot(border = TRUE, height = unit(2, "cm"), show_fraction = FALSE, ylim = NULL, axis = FALSE),
        MSI.status=(MS=MSI.MSI), Tissue.Type=(Type=MSI.Type), iCMS=(MS=MSI.iCMS), col=cols,
        simple_anno_size = unit(1.0, "cm"),height = unit(5, "cm"), annotation_name_gp=gpar(fontsize=25, "cm"),
        annotation_legend_param=list(
            iCMS=list(title_gp=gpar(fontsize=25), labels_gp=gpar(fontsize=25), grid_width = unit(1, "cm")),
            MSI.status=list(title_gp=gpar(fontsize=25), labels_gp=gpar(fontsize=25), grid_width = unit(1, "cm")),
            Tissue.Type=list(title_gp=gpar(fontsize=25), labels_gp=gpar(fontsize=25), grid_width = unit(1, "cm")))
		),
        right_annotation = rowAnnotation(row_barplot = anno_oncoprint_barplot(border = TRUE, width = unit(3, "cm"), show_fraction = FALSE, ,axis_param = list(side = "bottom", labels_rot = 90), axis = FALSE)),
		heatmap_legend_param = heatmap_legend_param, 
		height=unit(35,"cm"),
		width=unit(50,"cm"),
		pct_gp= grid::gpar(fontsize = 25, "cm") #Left percent font size
)
all.equal(colnames(MSI.crc), colnames(ht_list@matrix))

#pdf(paste0(work_dir,"01.Plot/pathway_sortMSI_230207.pdf"), width=30, height=20)
pdf(paste0(work_dir,"01.Plot/pathway_sortMSI_230303_tp53.pdf"), width=30, height=20)
print({draw(ht_list, merge_legend=TRUE, heatmap_legend_side = "left", annotation_legend_side = "bottom")})
dev.off()

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
p <- ggplot(myo, aes(x=RF.nearestCMS, y=Myofibroblasts)) + geom_boxplot(varwidth=TRUE)
q <- ggplot(rmna.myo, aes(x=RF.predictedCMS, y=Myofibroblasts)) + geom_boxplot()
#pdf(paste0(work_dir,"01.Plot/boxplot_nearest.pdf"), width=25, height=30)
pdf(paste0(work_dir,"01.Plot/boxplot_nearest.pdf"))
print({p})
dev.off()

pdf(paste0(work_dir,"01.Plot/boxplot_predicted.pdf"), width=25, height=30)
print({q})
dev.off()

#iCMS2, compare CMS2,CMS3 to CMS4
iCMS2_myo <- read.table("/data/keeyoung/gy_RNA/scCODA/scCODA_iCMS.txt", sep='\t', head=T)
iCMS2_myo <- iCMS2_myo[,c("Class","nearest","Myofibroblasts")]
iCMS2_myo <- iCMS2_myo[which(iCMS2_myo$Class %in% "iCMS2"),]
iCMS2_myo <- iCMS2_myo[-which(iCMS2_myo$nearest %in% "CMS1"),]
iCMS2_myo$nearest[which(iCMS2_myo$nearest %in% c("CMS2","CMS3"))]="others"
iCMS2_myo <- iCMS2_myo[,c("nearest","Myofibroblasts")]

piCMS2 <- ggplot(iCMS2_myo, aes(x=nearest, y=Myofibroblasts, color=nearest)) + ggtitle("iCMS2") + geom_boxplot(varwidth=TRUE, show.legend=FALSE) + 
	scale_color_manual(breaks=c("CMS4","others"), values=c("#31B404","#8080E5")) + 
#	theme(legend.position = "none") + 
	scale_y_continuous(limits = c(min(iCMS3_myo$Myofibroblasts), max(iCMS2_myo$Myofibroblasts))) +
	labs(x="CMS.Type") + theme_bw() + 
	theme(plot.title = element_text(hjust=0.5))
#pdf(paste0(work_dir,"01.Plot/boxplot_iCMS2_compare-CMS2,CMS3toCMS4.pdf"))
#print({piCMS2})
#dev.off()

#iCMS3, compare CMS2,CMS3 to CMS4
iCMS3_myo <- read.table("/data/keeyoung/gy_RNA/scCODA/scCODA_iCMS.txt", sep='\t', head=T)
iCMS3_myo <- iCMS3_myo[,c("Class","nearest","Myofibroblasts")]
iCMS3_myo <- iCMS3_myo[which(iCMS3_myo$Class %in% "iCMS3"),]
iCMS3_myo <- iCMS3_myo[-which(iCMS3_myo$nearest %in% "CMS1"),]
iCMS3_myo$nearest[which(iCMS3_myo$nearest %in% c("CMS2","CMS3"))]="others"
iCMS3_myo <- iCMS3_myo[,c("nearest","Myofibroblasts")]

piCMS3 <- ggplot(iCMS3_myo, aes(x=nearest, y=Myofibroblasts, color=nearest)) + ggtitle("iCMS3") + geom_boxplot(varwidth=TRUE) + 
	scale_color_manual(name="CMS.Type", breaks=c("CMS4","others"), values=c("#31B404","#8080E5")) + 
	scale_y_continuous(limits = c(min(iCMS3_myo$Myofibroblasts), max(iCMS2_myo$Myofibroblasts))) +
	labs(x="CMS.Type") + theme_bw() +
	theme(plot.title = element_text(hjust=0.5), axis.title.y=element_blank())

pdf(paste0(work_dir,"01.Plot/boxplot_iCMS3_compare-CMS2,CMS3toCMS4.pdf"))
print({piCMS3})
dev.off()

pdf(paste0(work_dir,"01.Plot/edit)boxplot_iCMS_compare-CMS2,CMS3toCMS4.pdf"))
#print({grid.arrange(piCMS2, piCMS3, shared_legend, ncol=2)})
print({grid.arrange(piCMS2, piCMS3,  ncol=2, widths=c(12,15))})
dev.off()


iCMS2_myo$nearest <- gsub("CMS4","iCMS2_CMS4", iCMS2_myo$nearest)
iCMS2_myo$nearest <- gsub("others","iCMS2_others", iCMS2_myo$nearest)

iCMS3_myo$nearest <- gsub("CMS4","iCMS3_CMS4", iCMS3_myo$nearest)
iCMS3_myo$nearest <- gsub("others","iCMS3_others", iCMS3_myo$nearest)

iCMS_myo <- rbind(iCMS2_myo, iCMS3_myo)
piCMS <- ggplot(iCMS_myo, aes(x=nearest, y=Myofibroblasts, color=nearest, fill=nearest)) + geom_boxplot(varwidth=TRUE)
pdf(paste0(work_dir,"01.Plot/boxplot_iCMS_compare-CMS2,CMS3toCMS4.pdf"))
#print({piCMS})
#print({piCMS + scale_color_manual(breaks=c("iCMS2_CMS4","iCMS2_others","iCMS3_CMS4","iCMS3_others"), values=c("#31B404","#8080E5","#31B404","#8080E5")) + 
#	scale_fill_manual(breaks=c("iCMS2_CMS4","iCMS2_others","iCMS3_CMS4","iCMS3_others"), values=c("#0040FF","#0040FF","#F5A9D0","#F5A9D0"))
#})
print({piCMS})
print({piCMS + scale_color_manual(breaks=c("iCMS2_CMS4","iCMS2_others","iCMS3_CMS4","iCMS3_others"), values=c("#0040FF","#0040FF","#F5A9D0","#F5A9D0")) + 
   scale_fill_manual(breaks=c("iCMS2_CMS4","iCMS2_others","iCMS3_CMS4","iCMS3_others"), values=c("#31B404","#8080E5","#31B404","#8080E5"))
})


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

