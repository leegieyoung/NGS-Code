#made by Giyoung Lee, If you want to help me,contact gy48085@gmail.com
source_dir="/data/keeyoung/scRNA/iCD/Code/source/"
source(paste0(source_dir,"source_sctype_cell_only_iCD.R"))

iCD.combined <- readRDS(paste0(pre_dir, "output/1.QC/","iCD.combined" , "_rds"))

#It reduced the difference between two groups
DefaultAssay(iCD.combined) <- "integrated"
QC.iCD.combined <- subset(iCD.combined, subset=nCount_RNA >5 & nCount_RNA < 2000)
#QC.iCD.combined <- iCD.combined
DefaultAssay(QC.iCD.combined) <- "integrated"
#Idents(QC.iCD.combined) <- QC.iCD.combined$integrated_snn_res.0.4
tissue = "iCD"
#for (i in seq(2.0,2.5,0.1)){
#	sctype_cell("iCD",QC.iCD.combined, i)
#}
sctype_cell("iCD",QC.iCD.combined, 0.8)
