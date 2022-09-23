#made by Giyoung Lee, If you want to help me,contact gy48085@gmail.com
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
source_dir="/data/keeyoung/scRNA/iCD/Code/source/"
source(paste0(source_dir,"source_barplot.R")
#source("/data/keeyoung/scRNA/iCD/Code/GIMAT_heatmap.R")
lapply(c("dplyr","Seurat","HGNChelper","openxlsx","scater","patchwork","metap","limma","multtest","ggplot2","cowplot","future.apply","ggraph","igraph","tidyverse","data.tree"), library, character.only = T)
pre_dir <- '/data/keeyoung/scRNA/iCD/'
work_dir <- paste0(pre_dir,"output/sctype/")
#db_ = paste0(pre_dir, "output/sctype/", "iCD_metadata.xlsx");
db_ = paste0(pre_dir, "output/sctype/", "iCD_metadata_220923.xlsx");
#GIMATs <- c("INF.macs","IgG","Activated DC","Activated Fibroblasts","ACKR1")
#function
sctype_Bcell <- function(Cell, Cell.combined, Change){
	work_dir <- paste0(pre_dir,"output/sctype/")
	if (!dir.exists(paste0(work_dir,Cell,"/"))){
	    dir.create(paste0(work_dir,Cell,"/"))
	}
	if (!dir.exists(paste0(work_dir,Cell,"/",Change))){
        dir.create(paste0(work_dir,Cell,"/",Change))
    }
	work_dir <- paste0(work_dir,Cell,"/",Change,"/")
	DefaultAssay(Cell.combined) <- "integrated"
	Cell.combined <- FindNeighbors(Cell.combined, dims = 1:13, verbose = FALSE)
	Cell.combined <- FindClusters(Cell.combined, resolution = Change, verbose = FALSE)
	Cell.combined <- RunUMAP(object = Cell.combined, dims = 1:13, verbose = FALSE)

	gs_list = gene_sets_prepare(db_, tissue)
	es.max = sctype_score(scRNAseqData = Cell.combined[["integrated"]]@scale.data, scaled = TRUE,
		                  gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

	cL_resutls = do.call("rbind", lapply(unique(Cell.combined@meta.data$seurat_clusters), function(cl){
	    es.max.cl = sort(rowSums(es.max[ ,rownames(Cell.combined@meta.data[Cell.combined@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
	    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(Cell.combined@meta.data$seurat_clusters==cl)), 10)
	}))
	sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

	# set low-confident (low ScType score) clusters to "unknown"
	sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
	print(sctype_scores[,1:3])

	Cell.combined@meta.data$customclassif = ""
	for(j in unique(sctype_scores$cluster)){
	  cl_type = sctype_scores[sctype_scores$cluster==j,];
	  Cell.combined@meta.data$customclassif[Cell.combined@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
	}

	png(filename=paste0(work_dir, "customclassif_res", Change, ".png"), width=200, height=54, units="cm", res=200)
	p1 <- (DimPlot(Cell.combined, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif',  label.size = 15, pt.size = 4, sizes.highlight = 1) +
	guides(color = guide_legend(override.aes = list(size=20), ncol=1,title.theme= element_text(size=20) )) + theme(legend.text = element_text(size=15)))

	p2 <- (DimPlot(Cell.combined, reduction = "umap", group.by = "group",label = TRUE, repel = TRUE, label.size = 30, pt.size = 4, sizes.highlight = 1) +
	guides(color = guide_legend(override.aes = list(size=20), ncol=1,title.theme= element_text(size=20) )) + theme(legend.text = element_text(size=15)))

	print({p1+p2})
	dev.off()

	cL_resutls=cL_resutls[order(cL_resutls$cluster),]; edges = cL_resutls; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL

	# prepare nodes
	nodes_lvl1 = sctype_scores[,c("cluster", "ncells")]; nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c();
	ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3",'#33332E','#A8A780','#5C6A0F','#3D5F16','#ff7f77','#fb9a00',"#33a99d","#15aa0f","#25aa0f","#35aa0f","#45aa0f","#55aa0f","#65aa0f","#12aa0f","#13aa0f","#14aa0f","#15aa0f","#16aa0f","#17aa0f","#18aa0f","#19aa0f","#20aa0f","#21aa0f","#22aa0f","#23aa0f","#24aa0f","#25aa0f","#26aa0f","#27aa0f","#28aa0f","#29aa0f","#30aa0f","#31aa0f","#32aa0f","#33aa0f","#34aa0f","#35aa0f","#36aa0f","#37aa0f","#38aa0f","#39aa0f","#40aa0f","#41aa0f","#42aa0f","#43aa0f","#44aa0f","#45aa0f","#46aa0f","#47aa0f","#48aa0f","#49aa0f","#50aa0f","#51aa0f","#52aa0f","#53aa0f","#54aa0f","#55aa0f","#56aa0f","#57aa0f","#58aa0f","#59aa0f")
	for (i in 1:length(unique(cL_resutls$cluster))){
	  dt_tmp = cL_resutls[cL_resutls$cluster == unique(cL_resutls$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
	}
	nodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
	files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
	nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]

	mygraph <- graph_from_data_frame(edges, vertices=nodes)

	gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) +
	  geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
	  theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="black", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1, label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")

	print("multiplot")
	png(filename=paste0(work_dir, "multiplot_res", Change, ".png"), width=150, height=54, units="cm", res=200)
	print({scater::multiplot(DimPlot(Cell.combined, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss, label.size = 20, pt.size = 4, sizes.highlight = 1), gggr, cols = 2)})
	DimPlot(Cell.combined, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss, label.size = 20, pt.size = 4, sizes.highlight = 1)
	print("Start DimPlot")
	DimPlot(Cell.combined, reduction = "umap", label = TRUE, repel = TRUE, label.size = 20, pt.size = 4, sizes.highlight = 1)
	dev.off()
	print("Barplot")
	png(filename=paste0(work_dir, "barplot_res", Change, ".png"), width=150, height=120, units="cm", res=200)
	print({Barplot(Cell,Cell.combined)})
	dev.off()
	
	#DoHeatmap
	#all.markers <- FindAllMarkers(Cell.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
	#all.markers %>% group_by(cluster) -> heatmap.allmarkers
	#DefaultAssay(Cell.combined) <- "RNA"
    #all.genes <- rownames(Cell.combined)
    #Cell.combined <- ScaleData(Cell.combined, verbose = FALSE, features=all.genes)
	#pdf(paste0(work_dir, "Cell_Doheatmap_",Change,"_RNA.pdf"), width=54, height=120)
	#print({
    #DoHeatmap(Cell.combined, features=heatmap.allmarkers$gene, size=20) + NoLegend() + theme(axis.text.y = element_text(size =12) )
    #})
	#dev.off()
	saveRDS(Cell.combined,paste0(work_dir,Cell,"_",Change,".rds"))
	#GIMAT_heatmap(Cell, Cell.combined, Change)
	return(Cell.combined)
}
