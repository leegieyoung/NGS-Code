#made by Giyoung Lee, If you want to help me,contact gy48085@gmail.com
lapply(c("dplyr","ggplot2"), library, character.only = T) 
MNP <- c("INF.macs"="#C4F875", "Activated DC"="#EEC053", "pDC"="#4FB0F2", "DC2"="#E37F67", "DC1"="#D5AFF9", "RES.macs"="#FDFAD1", "moDC"="#9EECF1")
Plasma <- c("IgG"="#E37F67","IgA"="#FDFAD1","IgM"="#9EECF1","Naive B"="#EEC053", "Memory B"="#4FB0F2")
#iCD <- c("T cells"="#4FB0F2","ILC"="#D5AFF9","Plasma cells"="#E37F67","B cells"="#EEC053","MNP"="#C4F875","pDC"="#FFFF00","Mast cells"="#E6E6E6","Stromal/glia"="#D7A651","cell cycling"="#0095FA","immunoregulation"="#148DDE","Naive/CM"="#1D78B5","CD8/cytotoxic"="#4686B0","resident memory"="#638eab")
Tcell <- c("Highly Activated T cells"="#D6D8DA","Tregs"="#FCD2FA","Naive T cells"="#C4F800","CM T cells"="#C4F875","Cytotoxic T cells"="#EEC053","Gamma delta T cell"="#D5AFF9","TFH-like"="#4FB0F2","CD8 Trm"="#FDFAD1","Remaining Trm"="#E37F67")
Bcell <- c("Nave B cells"="#D6D8DA","Memory B cells"="#FCD2FA")
Stromal <- c("CD36+EC"="#D5AFF9","ACKR1"="#C4F875","Pericytes"="#EEC053","Smooth Muscle"="#FDFAD1","Fibroblasts"="#E37F67","Activated Fibroblasts"="#FCD2FA","Glial cell"="#9EECF1")
iCD <- c("T cells"="#4FB0F2","ILC"="#D5AFF9","Plasma cells"="#E37F67","B cells"="#EEC053","MNP"="#C4F875","pDC"="#FFFF00","Mast cells"="#E6E6E6","Stromal/glia"="#D7A651")
iCD <- c("T cells"="#4FB0F2","Plasma cells"="#E37F67","B cells"="#EEC053","MNP"="#C4F875","pDC"="#F7F800","Mast cells"="#E6E6E6","Stromal/glia"="#D7A651","cell cycling"="#4c82a6","immunoregulation"="#1b4663","Naive/CM"="#0a304a","CD8/cytotoxic"="#054e80","resident memory"="#0d314a","ILC"="#502159")
MNP_barplot <- c("INF.macs","Activated DC","pDC","DC1","DC2","RES.macs","moDC")
Barplot <- function(Cell, samples){

if(Cell=="Plasma"){
Values=Plasma
} else if (Cell=="MNP"){
	Values=MNP
} else if (Cell=="iCD"){
	Values=iCD
} else if (Cell=="Tcell"){
	Values=Tcell
} else if (Cell=="Stromal"){
	Values=Stromal
} else if (Cell=="Bcell"){
    Values=Bcell
}

Infla <- samples[,samples@meta.data$group == "Infla"]
Uninf <- samples[,samples@meta.data$group == "Uninf"]
print("p1")
p1 <- samples@meta.data %>%
    dplyr::group_by(group, customclassif) %>%
    dplyr::count() %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(percent=100*n/sum(n)) %>%
    dplyr::ungroup() %>%
    ggplot(aes(x=group, y=percent, fill=customclassif)) + geom_col(colour= "black") + geom_bar(stat="identity", color="black") +
    theme(axis.title=element_text(size=120, face="bold"), axis.text.x = element_text(size=100, face = "bold"), axis.text.y = element_text(size =100, face = "bold"), legend.title = element_text(color = "black", size = 90, face = "bold"),
        legend.text = element_text(color = "black", size = 80)) +
#   scale_fill_brewer(palette = "Pastel1")
	scale_fill_manual(values=Values)

print("Infla")
p2 <- Infla@meta.data %>%
    dplyr::group_by(patients, customclassif) %>%
    dplyr::count() %>%
    dplyr::group_by(patients) %>%
    dplyr::mutate(percent=100*n/sum(n)) %>%
    dplyr::ungroup() %>%
    ggplot(aes(x=patients, y=percent, fill=customclassif)) + geom_col(colour= "black") + geom_bar(stat="identity", color="black") +
    theme(axis.title=element_text(size=120, face="bold"), axis.text.x = element_text(size=50, face = "bold",angle=90), axis.text.y = element_text(size =100, face = "bold"), legend.title = element_text(color = "black", size = 90, face = "bold"),
        legend.text = element_text(color = "black", size = 80)) +
#   scale_fill_brewer(palette = "Pastel1")
    scale_fill_manual(values=Values) +
	scale_x_discrete(limits = c("GSM3972009_69","GSM3972016_138","GSM3972013_128","GSM3972022_187","GSM3972020_181","GSM3972026_193","GSM3972028_196","GSM3972017_158","GSM3972024_190"))

print("Uninf")
p3 <- Uninf@meta.data %>%
    dplyr::group_by(patients, customclassif) %>%
    dplyr::count() %>%
    dplyr::group_by(patients) %>%
    dplyr::mutate(percent=100*n/sum(n)) %>%
    dplyr::ungroup() %>%
    ggplot(aes(x=patients, y=percent, fill=customclassif)) + geom_col(colour= "black") + geom_bar(stat="identity", color="black") +
    theme(axis.title=element_text(size=120, face="bold"), axis.text.x = element_text(size=50, face = "bold",angle=90), axis.text.y = element_text(size =100, face = "bold"), legend.title = element_text(color = "black", size = 90, face = "bold"),
        legend.text = element_text(color = "black", size = 80)) +
#   scale_fill_brewer(palette = "Pastel1")
    scale_fill_manual(values=Values) +
	scale_x_discrete(limits = c("GSM3972010_68","GSM3972015_135","GSM3972014_129","GSM3972021_186","GSM3972019_180","GSM3972025_192","GSM3972027_195","GSM3972018_159","GSM3972023_189")) 

p1/(p2+p3)
}
