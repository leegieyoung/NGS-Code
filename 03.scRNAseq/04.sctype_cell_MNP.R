#made by Giyoung Lee, If you want to help me,contact gy48085@gmail.com
source_dir="/data/keeyoung/scRNA/iCD/Code/source/"
source(paste0(source_dir,"source_sctype_cell.R"))
#iCD.combined <- readRDS(paste0(pre_dir, "output/1.QC/","iCD.combined" , "_rds"))

#It reduced the difference between two groups
#DefaultAssay(iCD.combined) <- "integrated"
#QC.iCD.combined <- subset(iCD.combined, subset=nCount_RNA >5 & nCount_RNA < 2000)
#DefaultAssay(QC.iCD.combined) <- "integrated"
#Idents(QC.iCD.combined) <- QC.iCD.combined$integrated_snn_res.0.8

#MNP
QC.iCD.combined <- readRDS(paste0(pre_dir, "output/sctype/iCDfor_iCD/0.8/","iCD.combined_0.8.rds"))
DefaultAssay(QC.iCD.combined) <- "integrated"
#Tcell
Idents(QC.iCD.combined) <- QC.iCD.combined$customclassif
MNP.combined <- subset(QC.iCD.combined, ident=c("MNP"))
DefaultAssay(MNP.combined) <- "integrated"
tissue = "MNP"

#for (i in seq(0.2,2.5,0.1)){
#	sctype_cell("MNP",MNP.combined, i)
#}

#sctype_cell("MNP", MNP.combined, 1.8)
sctype_cell("MNP", MNP.combined, 4.0)
