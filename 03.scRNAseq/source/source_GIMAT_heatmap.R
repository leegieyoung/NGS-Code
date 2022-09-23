#made by Giyoung Lee, If you want to help me,contact gy48085@gmail.com
lapply(c("dplyr","Seurat","HGNChelper","openxlsx","scater","patchwork","metap","limma","multtest","ggplot2","cowplot","future.apply","ggraph","igraph","tidyverse","data.tree"), library, character.only = T)
pre_dir <- '/data/keeyoung/scRNA/iCD/'
work_dir <- paste0(pre_dir,"output/sctype/")
#MNP_GIMAT <- c("INF.macs")
MNP_GIMAT <- c("Activated DC")
Tcell_GIMAT <- "Highly Activated T cells"
#Stromal_GIMAT <- c("Activated Fibroblasts")
Stromal_GIMAT <- c("ACKR1")
Plasma_GIMAT <- "IgG"

ADC <- c("Activated DC")
INF <- c("INF.macs")
AF <- c("Activated Fibroblasts")
ACKR1 <- c("ACKR1")

#Plasma_Vln <- c("IGHG1","IGHG2","IGHGP","IGHG3","IGHG4","MZB1","IGKC","JCHAIN","IGHA1","IGHA2","MTRNR2L12","MTRNR2L8","TMEM107","IGHM")
#MNP_Vln <- c("FTH1","C15orf48","SOD2","CD63","GAPDH","CCL3","CSTB","TXN","PFN1","HLA-DRA","CXCL8","IL1RN","G0S2","VIM","S100A11","HLA-DRB5","IDO1","CD74","HLA-DQA2","INHBA","NINJ1","FCER1G","LDHA","RGCC","BCL2A1","CXCL3","IGHA1","THBS1","JCHAIN")
#Stromal_Vln <- c("GAPDH","TMSB10","IFI27","IGHG3","DUSP23","IGHA1","JCHAIN","PTPRB","SOX6","IGHA1")
#Tcell_Vln <- c("HLA-DRA","CD38","IL2RA","CD25","CD40LG","GZMA","GZMB","PRF1")

IgG_Vln <- as.character(as.matrix(read.table(paste0(work_dir,"IgG_fdr_exist_expression_on_Cibersortx.txt"))))

ADC_Vln <- as.character(as.matrix(read.table(paste0(work_dir,"ADC_fdr_exist_expression_on_Cibersortx.txt"))))
INF_Vln <- as.character(as.matrix(read.table(paste0(work_dir,"INF_fdr_exist_expression_on_Cibersortx.txt"))))
HAT_Vln <- c("POU2AF1","IGHA1","XIST","JCHAIN")
ACKR1_Vln <- as.character(as.matrix(read.table(paste0(work_dir,"ACKR1_fdr_exist_expression_on_Cibersortx.txt"))))
AF_Vln <- as.character(as.matrix(read.table(paste0(work_dir,"AF_fdr_exist_expression_on_Cibersortx.txt"))))

Plasma_Vln <- IgG_Vln
Tcell_Vln <- HAT_Vln
#MNP_Vln <- INF_Vln
MNP_Vln <- ADC_Vln
#Stromal_Vln <- AF_Vln
Stromal_Vln <- ACKR1_Vln

#function

GIMAT_heatmap <- function(Cell, Cell.combined, Change){
	print("Make_dir")
    work_dir <- paste0(pre_dir,"output/sctype/")
    if (!dir.exists(paste0(work_dir,Cell,"/"))){
        dir.create(paste0(work_dir,Cell,"/"))
    }
    if (!dir.exists(paste0(work_dir,Cell,"/",Change))){
        dir.create(paste0(work_dir,Cell,"/",Change))
    }
	print("GIMAT modules")
	if(Cell=="Plasma"){
		Cell_GIMAT=Plasma_GIMAT
		Vln_GIMAT=Plasma_Vln
	} else if (Cell=="MNP"){
		Cell_GIMAT=MNP_GIMAT
		Vln_GIMAT=MNP_Vln
	} else if (Cell=="iCD"){
		Cell_GIMAT=iCD_GIMAT
	} else if (Cell=="Tcell"){
		Cell_GIMAT=Tcell_GIMAT
		Vln_GIMAT=Tcell_Vln
	} else if (Cell=="Stromal"){
		Cell_GIMAT=Stromal_GIMAT
		Vln_GIMAT=Stromal_Vln}

	if (!dir.exists(paste0(work_dir,Cell,"/",Change,"/",Cell_GIMAT))){
        dir.create(paste0(work_dir,Cell,"/",Change,"/",Cell_GIMAT))
    }
	Vln_dir <- paste0(work_dir,Cell,"/",Change,"/")
	work_dir <- paste0(work_dir,Cell,"/",Change,"/",Cell_GIMAT,"/")
	print({work_dir})
	print("Idents customclassif")
	Idents(Cell.combined) <- Cell.combined$customclassif
	#VlnPlot for Total samples
	DefaultAssay(Cell.combined) <- "RNA"
	png(filename=paste0(Vln_dir, "VlnPlot_",Cell_GIMAT,"_res", Change, ".png"), width=100, height=54, units="cm", res=200)
	print({
	VlnPlot(Cell.combined, features=Vln_GIMAT)
	})
	dev.off()

	Cell.combined <- subset(Cell.combined, idents=Cell_GIMAT)
	
	GIMAT.combined <- Cell.combined[,grep("(No|)GIMAT", Cell.combined$GIMATs)]
	Idents(GIMAT.combined) <- GIMAT.combined$GIMATs

	print("GIMAT_barcodes")
	p1 <- rownames(GIMAT.combined@meta.data[which(GIMAT.combined@meta.data$GIMATs=="GIMAT"),])

	write.table(p1, paste0(work_dir,"GIMAT_length",length(p1),".txt"))
	print("NoGIMAT_barcodes")
	p2 <- rownames(GIMAT.combined@meta.data[which(GIMAT.combined@meta.data$GIMATs=="NoGIMAT"),])

	write.table(p2, paste0(work_dir,"NoGIMAT_length",length(p2),".txt"))
	#all.markers <- FindAllMarkers(GIMAT.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
	#all.markers <- FindAllMarkers(GIMAT.combined, only.pos = TRUE, logfc.threshold = 0.25)
	#all.markers <- FindAllMarkers(GIMAT.combined, only.pos = TRUE)
	all.markers <- FindAllMarkers(GIMAT.combined, only.pos = FALSE)
	write.table(all.markers, paste0(work_dir,"all.markers_",Cell_GIMAT,".csv"),quote=F, sep=',', col.names=T)
	
	#DimPlot to compare
	png(filename=paste0(work_dir, Cell_GIMAT,"_res", Change, ".png"), width=100, height=54, units="cm", res=200)
	print({DimPlot(GIMAT.combined, reduction = "umap", group.by = "GIMATs",label = TRUE, repel = TRUE, label.size = 15, pt.size = 6, sizes.highlight = 1) +
		ggtitle(Cell_GIMAT) +
	    #guides(color = guide_legend(override.aes = list(size=20), ncol=1,title.theme= element_text(size=30) )) + theme(legend.text = element_text(size=30), plot.title=element_text(hjust=0.5, size=30))
		#Legend figure size
		guides(color = guide_legend(override.aes = list(size=40))) +
		theme(legend.text = element_text(size=50), plot.title=element_text(hjust=0.5, size=60, face="bold"))
	})
	dev.off()
	#DoHeatmap
	all.markers %>% group_by(cluster) -> heatmap.allmarkers
	DefaultAssay(GIMAT.combined) <- "RNA"
	all.genes <- rownames(GIMAT.combined)
	GIMAT.combined <- ScaleData(GIMAT.combined, verbose = FALSE, features=all.genes)
	pdf(paste0(work_dir, "Cell_Doheatmap_",Change,"_RNA.pdf"), width=54, height=120)
	print({
	#DoHeatmap(GIMAT.combined, features=heatmap.allmarkers$gene, size=20) + NoLegend() + theme(axis.text.y = element_text(size =60) ) 
	DoHeatmap(GIMAT.combined, features=Vln_GIMAT, size=20) + NoLegend() + theme(axis.text.y = element_text(size =60) )
	})
	dev.off()
	#Vln Plot in GIMATS
    DefaultAssay(GIMAT.combined) <- "RNA"
    png(filename=paste0(work_dir, "VlnPlot_",Cell_GIMAT,"_res", Change, ".png"), width=100, height=54, units="cm", res=200)
    print({
    VlnPlot(GIMAT.combined, features=Vln_GIMAT)
	saveRDS(GIMAT.combined,paste0(work_dir,Cell_GIMAT,"_",Change,".rds"))
    })   
    dev.off()


}
