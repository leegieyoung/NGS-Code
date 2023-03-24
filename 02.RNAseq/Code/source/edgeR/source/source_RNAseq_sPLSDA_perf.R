suppressMessages({
library(edgeR)
library(limma)
library(RColorBrewer)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)
library(repr)
library(statmod)
library(GO.db)
library(ComplexHeatmap)
library(mixOmics)
})

QC_sPLSDA <- function(name,iCD){
#DIR
print("Make_DIR")
if (!dir.exists(paste0(pre_dir,name,"/"))){
	dir.create(paste0(pre_dir,name,"/"))
}
work_dir <- paste0(pre_dir,name,"/")

if (!dir.exists(paste0(work_dir,"01.Plot"))){
	dir.create(paste0(work_dir,"01.Plot"))
}
if (!dir.exists(paste0(work_dir,"02.Table"))){
    dir.create(paste0(work_dir,"02.Table"))
}
if (!dir.exists(paste0(work_dir,"03.Pathway"))){
    dir.create(paste0(work_dir,"03.Pathway"))
}

#Loading data
print("Loading data")
gtf <- read.table("/data/keeyoung/gy_RNA/REFERENCE/gtf/anno_Homo_sapiens.GRCh38.100.gtf", sep=" ", head=T)
#go_category
go_cg <- read.table("/data/keeyoung/gy_RNA/REFERENCE/GOterm/go_categories.txt", sep=",", head =F)


#my.contrasts <- makeContrasts(INvsNO=Infla-Uninf, levels=sample.design)

#INvsNO <- c("Infla","Uninf")
#wINvsNO <- which(sample.group %in% INvsNO)

#col_INvsNO <- list(cINvsNO=c("Infla"="blue","Uninf"="yellow"))

###
print("End of loading raw data")
cnt <- iCD
print("Do edgeR and Annotation using org.Hs")
cnt.lcpm <- cpm(cnt, log=TRUE)
y <- DGEList(cnt, group=sample.group, genes=rownames(cnt.lcpm),)

idfound <- rownames(y) %in% mappedRkeys(org.Hs.egENSEMBL)
y <- y[idfound,]
#Ensembl add function and extract proten_coding.
m <- match(rownames(y), gtf$ensembl_id)
Func <- gtf$function.[m]
y$genes <- cbind(y$genes, Func)
g <- grep("protein_coding", y$genes$Func)
y <- y[g,]
egENTREZID <- toTable(org.Hs.egENSEMBL)
m <- match(rownames(y), egENTREZID$ensembl_id)
y$genes$entrezID <- egENTREZID$gene_id[m]

#Entrezid -> Gene Symbol
egSYMBOL <- toTable(org.Hs.egSYMBOL)
m <- match(y$genes$entrezID, egSYMBOL$gene_id)
y$genes$symbol <- egSYMBOL$symbol[m]

#Remove LOC gene
remove <- grep("^LOC", y$genes$symbol, value=T)
y <- y[which(!y$genes$symbol %in% remove),]

#Remove Duplicated Genesymbols
o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]
y_order <- y[order(y$genes$symbol, decreasing=TRUE),]
d <- duplicated(y_order$genes$symbol)
y_uniq <- y[!d,]

#edgeR filter
keep <- filterByExpr(y_uniq)
y_uniq <- y_uniq[keep, , keep.lib.sizes=FALSE]
y_uniq <- y_uniq[filterByExpr(y_uniq),]
y_uniq <- calcNormFactors(y_uniq)
y_uniq <- estimateDisp(y_uniq,sample.design)
print({"EdgeR filter"})

#CPM
y_uniq.lcpm <- cpm(y_uniq, log=TRUE)
#write.table(y_uniq.lcpm, paste0(work_dir,"02.Table/cpm.csv"), sep='\t', col.names=T, row.names=T, quote=F)

#sPLSDA
y_uniq.lcpm.t <- t(y_uniq.lcpm)
y_uniq_lcpm.splsda <- splsda(y_uniq.lcpm.t, sample.group, ncomp=10)

#png(filename=paste0(work_dir,"01.Plot/sPLSDA_pca.png"), width=60, height=50, units="cm", res=200)
#print({plotIndiv(y_uniq_lcpm.splsda, comp = 1:2,
#        group = sample.group, ind.names = FALSE,
#        ellipse = TRUE,
#        legend = TRUE,
#        title="PLS-DA"
#)})
#dev.off()
##Tuning sPLSDA
#y_uniq_lcpm.perf.splsda <- perf(y_uniq_lcpm.splsda, validation = "Mfold", 
#                          folds = 10, nrepeat = 50, # use repeated cross-validation
#                          progressBar = TRUE, auc = TRUE, # include AUC values
#						  cpus = 36)

print({"Start readRDS"})
print({work_dir})
y_uniq_lcpm.perf.splsda <- readRDS(paste0(work_dir,"y_uniq_lcpm.perf.splsda.rds"))
print({"End of readRDS"})

print({"Start of make plot of perf"})
png(filename=paste0(work_dir,"01.Plot/sPLSDA_perf.png"), width=60, height=50, units="cm", res=200)
print({plot(y_uniq_lcpm.perf.splsda, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")})
dev.off()
print({"End of make plot of perf"})
#saveRDS(y_uniq_lcpm.perf.splsda,paste0(work_dir,"y_uniq_lcpm.perf.splsda.rds"))

#Selecting the number of variables
list.keepX <- c(1:10,  seq(20, 300, 10))
y_uniq_lcpm.tune.splsda <- tune.splsda(y_uniq.lcpm.t, sample.group, ncomp=3,
						validation = 'Mfold',
						folds= 10, nrepeat = 50,
						dist = 'max.dist',
						measure = 'overall',
						test.keepX = list.keepX,
						cpus = 36)
print({"Start of saveRDS"})
saveRDS(y_uniq_lcpm.tune.splsda,paste0(work_dir,"y_uniq_lcpm.tune.splsda.rds"))
#y_uniq_lcpm.tune.splsda <- readRDS(paste0(work_dir,"y_uniq_lcpm.tune.splsda.rds"))
print({"End of saveRDS"})
print({"Make plot of tune"})
png(filename=paste0(work_dir,"01.Plot/sPLSDA_tune.png"), width=60, height=50, units="cm", res=200)
print({plot(y_uniq_lcpm.tune.splsda, col = color.jet(3),
		size.xlabel=rel(2), size.ylabel=rel(2)
)})
dev.off()
print({"End of plot of tune"})

optimal.ncomp <- y_uniq_lcpm.tune.splsda$choice.ncomp$ncomp
#optimal.ncomp <- 3 
optimal.keepX <- y_uniq_lcpm.tune.splsda$choice.keepX[1:optimal.ncomp]

opti <- cbind(optimal.ncomp, optimal.keepX)
write.table(opti, paste0(work_dir,"optimal.keepX.txt"),sep='\t', quote=F)

y_uniq.final.splsda <- splsda(y_uniq.lcpm.t, sample.group,
						ncomp = optimal.ncomp,
						keepX = optimal.keepX)

#png(filename=paste0(work_dir,"01.Plot/sPLSDA_final_comp1-2.png"), width=60, height=50, units="cm", res=200)
#print({plotIndiv(y_uniq.final.splsda, comp = c(1,2),
#print({plotIndiv(y_uniq.final.splsda, comp=1,
#        group = sample.group, ind.names = FALSE,
#		ellipse = TRUE,
#		legend = TRUE,
#		title=paste0("sPLS-DA, comp ","& 1,2"),
#		size.title=rel(3),
#		size.xlabel=rel(2), size.ylabel=rel(2),
#		size.axis=rel(0.8),
#		size.legend.title=rel(2.3), size.legend=rel(2)
#		)})
#dev.off()

print({"Make plot of heatmap"})
sample.col <- sample.group
sample.col <- gsub("other","#F68B33",sample.col) #orange
sample.col <- gsub("NSNP","#388ECC",sample.col) #blue

#legend=list(legend=levels(sample.group),
legend=list(legend=c("other","NSNP"),
			col=c("#F68B33","#388ECC"),
			title="Last behavior",
			cex=1.4)

png(filename=paste0(work_dir,"01.Plot/sPLSDA_final_heatmap_comp1.png"), width=70, height=30, units="cm", res=200)
print({cim(y_uniq.final.splsda, row.sideColors = sample.col,
			keysize=c(0.2,0.9),
			legend = legend)
	})
dev.off()


}

