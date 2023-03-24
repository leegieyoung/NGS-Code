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
name="test"
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

Cell_frac <- read.table("/data/keeyoung/scRNA/cibersortx/output/02.Imputed_Cell_Fraction/IBD/IBD_cpm_221102/CIBERSORTx_Adjusted.txt", sep='\t')
colnames(Cell_frac) <- Cell_frac[1,]
Cell_frac <- Cell_frac[-1,]
rownames(Cell_frac) <- Cell_frac[,1]
Cell_frac <- Cell_frac[,-1]
Cell_frac <- Cell_frac[,c(-24, -25,-26)]
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

#Unknown <- grep("Unknown",sample.pheno)
#sample.pheno <- sample.pheno[-Unknown]

#sample.group <- factor(sample.pheno, levels=c("NSNP","stricture","penetrating"), order = T)
sample.group <- factor(sample.pheno)

MGI_CD_infla <- Cell_frac[match(names(sample.pheno), rownames(Cell_frac)),]
m_MGI_CD_infla <- as.matrix((MGI_CD_infla))
mode(m_MGI_CD_infla) <- "numeric"

#Heatmap
sc_MGI_CD_infla <- m_MGI_CD_infla
for(i in 1:length(rownames(sc_MGI_CD_infla))){
sc_MGI_CD_infla[i,] <- scale(sc_MGI_CD_infla[i,], center=TRUE, scale=FALSE)
}

pdf(paste0(work_dir,"01.Plot/Cell_frac_Heatmap.pdf"), width=25, height=30)
anno <- HeatmapAnnotation(iCMS=(iCMS=sample.pheno), CMS_frac=anno_barplot(cbind(sample.CMS2,sample.CMS3,sample.CMS1,sample.CMS4), gp=gpar(fill=c("blue","pink","yellow","green"))),
    MS=(MS=sample.MSI), Type=(Type=sample.Type), col=cols,
    simple_anno_size = unit(2.5, "cm"),height = unit(20, "cm"), annotation_name_gp=gpar(fontsize=35, fontface="bold"),
    annotation_legend_param=list(iCMS=list(title_gp=gpar(fontsize=35, fontface="bold"), labels_gp=gpar(fontsize=30)), MS=list(title_gp=gpar(fontsize=35, fontface="bold"), labels_gp=gpar(fontsize=30)), Type=list(title_gp=gpar(fontsize=35, fontface="bold"), labels_gp=gpar(fontsize=30)))
)

print({Heatmap(y_uniq_vali.lcpm.scale, top_annotation=anno,
        column_names_gp=grid::gpar(fontsize=25), row_names_gp = grid::gpar(fontsize = 0.1),
        row_order = Total_genes,
        width=unit(42, "cm"), height=unit(25,"cm"),
        heatmap_legend_param = list(title = "Mat", labels_gp = gpar(fontsize = 30), title_gp = gpar(fontsize = 30, fontface="bold"), border="red",grid_width = unit(1.5, "cm"))
        #column_title = "715_ordered genes(used 674 genes)", column_title_gp = gpar(fontsize = 45, fontface = "bold")
)})
dev.off()

#UnCentered
UnCenteredEuclidean <- function(x,y) {
        1 - ( x/sqrt( x %*% x ) ) %*% ( y/sqrt( y %*% y ) )
}
pdf(paste0(work_dir,"01.Plot/UnCentered_Cell_frac_Heatmap.pdf"), width=25, height=30)
tm_MGI_CD_infla <- t(m_MGI_CD_infla)
Heatmap(tm_MGI_CD_infla, clustering_distance_rows = UnCenteredEuclidean, clustering_distance_column = UnCenteredEuclidean,
	column_names_gp=grid::gpar(fontsize=25), row_names_gp = grid::gpar(fontsize = 25),
	column_order=names(sort(sample.pheno)),
	width=unit(42, "cm"), height=unit(42,"cm") )
dev.off()

source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_normalized_ensemblToGenesymbol_CD_in_Article_246genes.R")) ##221130
RNAseq("CD_in_Total_article",iCD)
