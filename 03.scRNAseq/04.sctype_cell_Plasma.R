#made by Giyoung Lee, If you want to help me,contact gy48085@gmail.com
source_dir="/data/keeyoung/scRNA/iCD/Code/source/"
source(paste0(source_dir,"source_sctype_cell.R")

#iCD.combined <- readRDS(paste0(pre_dir, "output/1.QC/","iCD.combined" , "_rds"))

#It reduced the difference between two groups
#DefaultAssay(iCD.combined) <- "integrated"
#QC.iCD.combined <- subset(iCD.combined, subset=nCount_RNA >5 & nCount_RNA < 2000)
#DefaultAssay(QC.iCD.combined) <- "integrated"
#Idents(QC.iCD.combined) <- QC.iCD.combined$integrated_snn_res.0.8

#
QC.iCD.combined <- readRDS(paste0(pre_dir, "output/sctype/iCDfor_iCD/0.8/","iCD.combined_0.8.rds"))
DefaultAssay(QC.iCD.combined) <- "integrated"
#Tcell
Idents(QC.iCD.combined) <- QC.iCD.combined$customclassif

#Plasma
Plasma.combined <- subset(QC.iCD.combined, ident=c("Plasma cells"))
DefaultAssay(Plasma.combined) <- "integrated"
tissue = "Plasma"

#for (i in seq(0.2,1.5,0.1)){
#	sctype_cell("Plasma", Plasma.combined, i)
#}

#for (i in seq(1.5,4.0,0.5)){
#	sctype_cell("Plasma", Plasma.combined, i)
#}


sctype_cell("Plasma", Plasma.combined, 4.0)
sctype_cell("Plasma", Plasma.combined, 6.0)
