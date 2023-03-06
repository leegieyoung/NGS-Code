setwd("/data/keeyoung/scRNA/iCD/output/sctype")
QC.iCD.combined$celltype <- ""
QC.iCD.combined$celltype <- gsub("cell cycling","Tcell",QC.iCD.combined$celltype)
QC.iCD.combined$celltype <- gsub("immunoregulation","Tcell",QC.iCD.combined$celltype)
QC.iCD.combined$celltype <- gsub("CD8/cytotoxic","Tcell",QC.iCD.combined$celltype)
QC.iCD.combined$celltype <- gsub("resident memory","Tcell",QC.iCD.combined$celltype)
QC.iCD.combined$celltype <- gsub("Naive/CM","Tcell",QC.iCD.combined$celltype)

p4 <- QC.iCD.combined@meta.data %>% dplyr::group_by(patients, group,celltype) %>% dplyr::count() %>% dplyr::group_by(patients) %>% dplyr::mutate(percent=100*n/sum(n)) %>% dplyr::ungroup()

p4$patients <- gsub("GSM3972013_128","pat.7.inf",p4$patients)
p4$patients <- gsub("GSM3972016_138","pat.8.inf",p4$patients)
p4$patients <- gsub("GSM3972020_181","pat.11.inf",p4$patients)
p4$patients <- gsub("GSM3972022_187","pat.12.inf",p4$patients)
p4$patients <- gsub("GSM3972017_158","pat.10.inf",p4$patients)
p4$patients <- gsub("GSM3972024_190","pat.13.inf",p4$patients)
p4$patients <- gsub("GSM3972026_193","pat.14.inf",p4$patients)
p4$patients <- gsub("GSM3972028_196","pat.15.inf",p4$patients)
p4$patients <- gsub("GSM3972014_129","pat.7.uninf",p4$patients)
p4$patients <- gsub("GSM3972015_135","pat.8.uninf",p4$patients)
p4$patients <- gsub("GSM3972019_180","pat.11.uninf",p4$patients)
p4$patients <- gsub("GSM3972021_186","pat.12.uninf",p4$patients)
p4$patients <- gsub("GSM3972018_159","pat.10.uninf",p4$patients)
p4$patients <- gsub("GSM3972023_189","pat.13.uninf",p4$patients)
p4$patients <- gsub("GSM3972025_192","pat.14.uninf",p4$patients)
p4$patients <- gsub("GSM3972027_195","pat.15.uninf",p4$patients)
p4$patients <- gsub("GSM3972009_69","pat.5.inf",p4$patients)
p4$patients <- gsub("GSM3972010_68","pat.5.uninf",p4$patients)

p5 <-p4[order(p4$patients, decreasing=F),]

Values <- c("Tcell"="#0040FF","Mast cells"="#BDBDBD","Stromal/glia"="#8A4B08","Plasma cells"="#FF0000","B cells"="#FE9A2E","MNP"="#3ADF00")

png("celltype_patients.png",  width=20, height=12, units="cm", res=200)
p6 <- p5 %>% ggplot(aes(x=patients, y=percent, fill=celltype)) + geom_col(colour= "black") + geom_bar(stat="identity", color="black") + theme(axis.text.x=element_text(angle=90, hjust=1)) + scale_fill_manual(values=Values)
print({p6})
dev.off()

png("celltype_group.png",  width=20, height=12, units="cm", res=200)
p1 <- QC.iCD.combined@meta.data %>% dplyr::group_by(group,celltype) %>% dplyr::count() %>% dplyr::group_by(group) %>% dplyr::mutate(percent=100*n/sum(n)) %>% dplyr::ungroup() %>% ggplot(aes(x=group, y=percent, fill=celltype)) + geom_col(colour= "black") + geom_bar(stat="identity", color="black") + theme(axis.text.x=element_text(angle=90, hjust=1)) + scale_fill_manual(values=Values)
print({p1})
dev.off()


#
Idents(QC.iCD.combined) <- QC.iCD.combined$celltype
subtype <- QC.iCD.combined$celltype
Gene <- rownames(GetAssayData(QC.iCD.combined, slot='data',assay='integrated'))
QC.iCD.counts <- GetAssayData(QC.iCD.combined, slot='counts',assay='RNA')[Gene,]
QC.iCD.counts <- as.matrix(x=QC.iCD.counts)

colnames(QC.iCD.counts) <- subtype

write.table(QC.iCD.counts, "celltype_group.txt", sep='\t', col.names=T)
