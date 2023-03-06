#made by Giyoung Lee, If you want to help me,contact gy48085@gmail.com
library(Seurat)
library(patchwork)
library(metap)
library(limma)
library(multtest)
library(ggplot2)
library(cowplot)
library(future.apply)
library(dplyr)

pre_dir <- '/data/keeyoung/scRNA/iCD/'

#--------------integeration
Sample_list <- c("GSM3972009_69","GSM3972013_128","GSM3972016_138","GSM3972020_181","GSM3972022_187","GSM3972017_158","GSM3972024_190","GSM3972026_193","GSM3972028_196", "GSM3972010_68", "GSM3972014_129","GSM3972015_135", "GSM3972019_180", "GSM3972021_186", "GSM3972018_159", "GSM3972023_189", "GSM3972025_192", "GSM3972027_195") 

Assign_list <- c("Infla5", "Infla7", "Infla8", "Infla11", "Infla12", "Infla10", "Infla13", "Infla14", "Infla15", "Uninf5", "Uninf7", "Uninf8", "Uninf11", "Uninf12", "Uninf10", "Uninf13", "Uninf14", "Uninf15")
for (i in 1:length(Assign_list)) {
	assign(Assign_list[i], readRDS(paste0(pre_dir,"output/1.QC/", Sample_list[i], '_rds')))
}

#-----------------merge
cell.ids <- c("C69","C128","C138","C181","C187","C158", "C190", "C193", "C196", "C68", "C129", "C135", "C180", "C186", "C159", "C189", "C192", "C195") 
#iCD <- merge(Infla5, y=c(Infla7,Infla8, Infla11, Infla12, Infla10, Infla13, Infla14, Infla15, Uninf5, Uninf7, Uninf8, Uninf11, Uninf12, Uninf10, Uninf13, Uninf14, Uninf15), add.cell.ids = cell.ids, project = "iCD", merge.data = TRUE) 
iCD <- merge(eval(parse(text=Assign_list[1])), y=c(eval(parse(text = Assign_list[2])),
    eval(parse(text = Assign_list[3])),
    eval(parse(text = Assign_list[4])),
    eval(parse(text = Assign_list[5])),
    eval(parse(text = Assign_list[6])),
    eval(parse(text = Assign_list[7])),
    eval(parse(text = Assign_list[8])),
    eval(parse(text = Assign_list[9])),
    eval(parse(text = Assign_list[10])),
    eval(parse(text = Assign_list[11])),
    eval(parse(text = Assign_list[12])),
    eval(parse(text = Assign_list[13])),
    eval(parse(text = Assign_list[14])),
    eval(parse(text = Assign_list[15])),
    eval(parse(text = Assign_list[16])),
    eval(parse(text = Assign_list[17])),
    eval(parse(text = Assign_list[18]))),
    add.cell.ids = cell.ids, project = "iCD", merge.data = TRUE)

#-----------------integration
groups <- iCD@meta.data$orig.ident
names(groups) <- colnames(iCD)
iCD <- AddMetaData(object=iCD, metadata=groups, col.name="group")
iCD.list <- SplitObject(iCD, split.by ="group")

# normalize and identify variable features for each dataset independently
iCD.list <- lapply(X = iCD.list, FUN = function(x) {
	x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = iCD.list)

#It takes a longtime...
#Default l2.nrom =T , reduction  cca, normalization.method = LogRnormalize, featrues=2,000
iCD.anchors <- FindIntegrationAnchors(object.list = iCD.list, normalization.method = "LogNormalize", reduction = "cca", anchor.features = features)

iCD.combined <- IntegrateData(anchorset = iCD.anchors, dims=1:50)

DefaultAssay(iCD.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
iCD.combined <- ScaleData(iCD.combined, verbose = FALSE)
iCD.combined <- RunPCA(iCD.combined, verbose = FALSE)
iCD.combined <- JackStraw(iCD.combined, num.replicate = 100, dims=50)
iCD.combined <- ScoreJackStraw(iCD.combined, dims = 1:50)

#https://hbctraining.github.io/scRNA-seq/lessons/elbow_plot_metric.html
pct <- iCD.combined[["pca"]]@stdev / sum(iCD.combined[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co1

co2 <- sort(which((pct[1:length(pct) -1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2

pcs <- min(co1,co2)
pcs

plot_df <- data.frame(pct=pct, cumu = cumu, rank = 1:length(pct))

pdf("/data/keeyoung/scRNA/iCD/output/1.QC/ggplot_elbowplot.pdf")
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()

#Overlab Umap
pdf("/data/keeyoung/scRNA/iCD/output/1.QC/elbowplot.pdf")
ElbowPlot(iCD.combined, ndims=50)
dev.off()

pdf("/data/keeyoung/scRNA/iCD/output/1.QC/JackStrawPlot.pdf")
JackStrawPlot(iCD.combined, dims = 1:50)
dev.off()

saveRDS(iCD.combined, file = paste0(pre_dir, "output/1.QC/","raw.iCD.combined" , "_rds"))

