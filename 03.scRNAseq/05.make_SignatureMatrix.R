lapply(c("dplyr","Seurat","HGNChelper","openxlsx","scater","patchwork","metap","limma","multtest","ggplot2","cowplot","future.apply","ggraph","igraph","tidyverse","data.tree"), library, character.only = T)
pre_dir <- '/data/keeyoung/scRNA/iCD/'
work_dir <- paste0(pre_dir,"output/sctype/")

Change="4"

Plasma <- paste0(work_dir,"Plasma/",Change,"/","Plasma_",Change,".rds")
Tcell <- paste0("Tcell_",Change,".rds")
MNP <- paste0("MNP_",Change,".rds")
Stromal <- paste0("Stromal_",Change,".rds")

Cell_list <- c("Plasma","Tcell","MNP","Stromal","Bcell")
#Cell_list="Plasma"
Change="4"

merge=""
for (i in 1:length(Cell_list)){
	assign(Cell_list[i], readRDS(paste0(work_dir,Cell_list[i],"/",Change,"/",Cell_list[i],"_",Change,".rds")))
	Cell <- eval(parse(text=Cell_list[i]))
	Idents(Cell) <- Cell$customclassif
	subtype = Cell$customclassif
	#2000 Genes on integrated.
	#Gene <- rownames(GetAssayData(object=Cell, slot="data", assay="integrated"))
	#30000 over Genes on RNA.
	Gene <- rownames(GetAssayData(object=Cell, slot="data", assay="RNA"))
	Cell <- GetAssayData(object=Cell, slot="counts",assay="RNA")[Gene,]
	Cell <- as.matrix(x=Cell)
	print("Is it equal that barcode and Cell name ?")
	print(all.equal(colnames(Cell), names(subtype)))
	colnames(Cell) = subtype
	print(dim(Cell))
	merge <- cbind(merge, Cell)
}

All_Celltype_counts <- merge[,order(colnames(merge))]
write.table(All_Celltype_counts,paste0(work_dir,"All_Celltype_counts_220923_addBcell.txt"), sep="\t", col.names=T)

for (i in 1:length(Cell_list)){
    assign(Cell_list[i], readRDS(paste0(work_dir,Cell_list[i],"/",Change,"/",Cell_list[i],"_",Change,".rds")))
    Cell <- eval(parse(text=Cell_list[i]))
	DefaultAssay(Cell) <- 'integrated'
	Cell <- Cell[,grep("(No|)GIMAT", Cell$GIMATs)]
    Idents(Cell) <- Cell$customclassif
    #Merge <- merge(Merge, y=Cell, project= "iCD", merge.data=TRUE)
	assign(Cell_list[i], Cell)
}

p1 <- Merge@meta.data %>% dplyr::group_by(customclassif) %>% dplyr::count() %>% dplyr::group_by(customclassif) %>% dplyr::mutate(percent=100*n/sum(n))

Merge <- merge(Plasma, c(Tcell,MNP,Stromal), project= "iCD", merge.data=TRUE)
#saveRDS(Merge, file = paste0(work_dir,"Inflamed_iCD" , ".rds"))



