source_dir="/data/keeyoung/scRNA/iCD/Code/source/"
source(paste0(source_dir,"source_sctype_cell.R"))
Tcells <- c("immunoregulation","cell cycling","CD8/cytotoxic","resident memory","Naive/CM")
MNP <- "MNP"
Plasma <- "Plasma cells"
Stromal <- "Stromal/glia"
Bcells <- "B cells"
Infla <- c(
"GSM3972009_69",
"GSM3972013_128",
"GSM3972016_138",
"GSM3972020_181",
"GSM3972022_187",
"GSM3972017_158",
"GSM3972024_190",
"GSM3972026_193",
"GSM3972028_196"
)

Uninf <-c(
"GSM3972010_68",
"GSM3972014_129",
"GSM3972015_135",
"GSM3972019_180",
"GSM3972021_186",
"GSM3972018_159",
"GSM3972023_189",
"GSM3972025_192",
"GSM3972027_195"
)
#ConvertFC(Tcells,Tcells,QC.iCD.combined)
ConvertFC <- function(name,Cell, raw.combined){

Cell.combined <- subset(raw.combined, ident=Cell)
unique(Idents(Cell.combined))
Gene <- rownames(GetAssayData(object=Cell.combined, slot="data", assay="RNA"))

Idents(Cell.combined) <- Cell.combined$patients
Infla_merge=""
for ( i in 1:length(Infla)) {
print(Infla[i])
Infla.combined <- subset(Cell.combined, ident=Infla[i])
#print(head(Idents(Infla.combined)))
Infla.counts <- GetAssayData(Infla.combined, slot="counts", assay="RNA")[Gene,]
Infla.counts <- as.matrix(x=Infla.counts)
Infla.counts <- rowSums(Infla.counts)
Infla_merge <- cbind(Infla_merge, Infla.counts)
}
Infla_merge <- Infla_merge[,-1]

colnames(Infla_merge) <- Infla
print(head(Infla_merge))

Uninf_merge=""
for ( i in 1:length(Uninf)) {
print(Uninf[i])
Uninf.combined <- subset(Cell.combined, ident=Uninf[i])
#print(head(Idents(Uninf.combined)))
Uninf.counts <- GetAssayData(Uninf.combined, slot="counts", assay="RNA")[Gene,]
Uninf.counts <- as.matrix(x=Uninf.counts)
Uninf.counts <- rowSums(Uninf.counts)
Uninf_merge <- cbind(Uninf_merge, Uninf.counts)
}
Uninf_merge <- Uninf_merge[,-1]

colnames(Uninf_merge) <- Uninf
print(head(Uninf_merge))

result <- cbind(Infla_merge, Uninf_merge)
print(dim(result))
print(head(result))

write.table(result, paste0("/data/keeyoung/scRNA/iCD/output/convertFeatureCounts/",name,"_9patients.txt"), sep=' ', quote=F, col.names=T, row.names=T)
}

