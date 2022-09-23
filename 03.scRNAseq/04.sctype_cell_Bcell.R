#made by Giyoung Lee, If you want to help me,contact gy48085@gmail.com
source_dir="/data/keeyoung/scRNA/iCD/Code/source/"
source(paste0(source_dir,"source_sctype_cell_only_Bcell.R"))
#iCD.combined <- readRDS(paste0(pre_dir, "output/1.QC/","iCD.combined" , "_rds"))

#It reduced the difference between two groups
#DefaultAssay(iCD.combined) <- "integrated"
#QC.iCD.combined <- subset(iCD.combined, subset=nCount_RNA >5 & nCount_RNA < 2000)
#DefaultAssay(QC.iCD.combined) <- "integrated"
#Idents(QC.iCD.combined) <- QC.iCD.combined$integrated_snn_res.0.8

#
QC.iCD.combined <- readRDS(paste0(pre_dir, "output/sctype/iCDfor_iCD/0.8/","iCD.combined_0.8.rds"))
DefaultAssay(QC.iCD.combined) <- "integrated"
#Bcell
Idents(QC.iCD.combined) <- QC.iCD.combined$customclassif
Bcell.combined <- subset(QC.iCD.combined, ident=c("B cells"))
DefaultAssay(Bcell.combined) <- "integrated"
tissue = "Bcell"

#Only_Bcell
#for (i in seq(0.8,1.6,0.1)){
#	sctype_Bcell("Bcell", Bcell.combined, i)
#}

#for (i in seq(2.0,4.0,0.5)){
#	sctype_Bcell("Bcell", Bcell.combined, i)
#}

#sctype_Bcell("Bcell", Bcell.combined, 1.4)
sctype_Bcell("Bcell", Bcell.combined, 4.0)
