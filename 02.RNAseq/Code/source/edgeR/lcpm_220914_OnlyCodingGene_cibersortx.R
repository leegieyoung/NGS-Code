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

sample.meta <- read.table("/data/keeyoung/scRNA/iCD/01.Samplelist/sample.csv", sep=',')
#Extract
v2 <- grep("CD_infla",sample.meta$V2)
sample.meta <- sample.meta[v2,]
v4 <- grep("Totalseq", sample.meta$V4)
sample.meta <- sample.meta[v4,]
sample.ID <- sample.meta$V1
sample.ID <- gsub("-",".", sample.ID)
head(sample.meta)
sample.pheno <- sample.meta$V3
sample.pheno <- gsub("^R$","Likely NDR", sample.pheno)
sample.pheno <- gsub("^L$","Likely DR", sample.pheno)
names(sample.pheno) <- sample.ID
sample.group <- factor(sample.pheno)

y_uniq <- read.table("/data/keeyoung/scRNA/cibersortx/output/02.Imputed_Cell_Fraction/CD_infla/lcpm_220914_OnlyCodingGene/CIBERSORTx_Adjusted.txt", sep='\t')

#lcpm
y_uniq <- y_uniq[,c(-6,-19,-23,-24,-25)]
y_uniq_lcpm.splsda <- splsda(y_uniq, sample.group, ncomp=5)

work_dir <- '/data/keeyoung/gy_RNA/06.output/CD_in/04.sPLSDA/lcpm_220914_OnlyCodingGene/'
#PLSDA
png(filename=paste0(work_dir,"sPLSDA_plsda.png"), width=60, height=50, units="cm", res=200)
print({plotIndiv(y_uniq_lcpm.splsda, comp = 1:2,
        group = sample.group, ind.names = FALSE,
        ellipse = TRUE,
        legend = TRUE,
        title="PLS-DA",
        size.title=rel(3),
        size.xlabel=rel(2), size.ylabel=rel(2),
        size.axis=rel(0.8),
        size.legend.title=rel(2.3), size.legend=rel(2)
	
)})
dev.off()

#Perf
y_uniq_lcpm.perf.splsda <- perf(y_uniq_lcpm.splsda, validation = "Mfold", #Mfold, but mininal smaples should use LOO validation.
                          folds = 5, nrepeat = 100, # use repeated cross-validation
#						  folds=5, 
                          progressBar = TRUE, auc = TRUE) # include AUC values
saveRDS(y_uniq_lcpm.perf.splsda, paste0(work_dir,"sPLSDA_perf.rds"))

png(filename=paste0(work_dir,"sPLSDA_perf.png"), width=60, height=50, units="cm", res=200)
print({plot(y_uniq_lcpm.perf.splsda, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")})
dev.off()

#Check ncomp
y_uniq_lcpm.perf.splsda$choice.ncomp
#ncomp=1

list.keepX <- c(1:18)
y_uniq_lcpm.tune.splsda <- tune.splsda(y_uniq, sample.group, ncomp=1,
                        validation = 'Mfold',
#                        validation = 'loo',
                        folds= 5, nrepeat = 100,
                       dist = 'max.dist',
                       measure = 'overall',
                        test.keepX = list.keepX,
                        cpus = 36)
png(filename=paste0(work_dir,"sPLSDA_tune.png"), width=60, height=50, units="cm", res=200)
print({plot(y_uniq_lcpm.tune.splsda, col = color.jet(1),
        size.xlabel=rel(2), size.ylabel=rel(2)
)})
dev.off()

optimal.ncomp <- 1
#optimal.ncomp <- y_uniq_lcpm.tune.splsda$choice.ncomp$ncomp
optimal.keepX <- y_uniq_lcpm.tune.splsda$choice.keepX
#optimal.keepX <- y_uniq_lcpm.tune.splsda$choice.keepX[1]
#optimal.keepX <- c(1:5)
optimal <- cbind(optimal.ncomp, optimal.keepX)
write.table(optimal,paste0(work_dir,"keepX.txt"), sep='\t')
#sPLSDA
y_uniq.final.splsda <- splsda(y_uniq, sample.group,
                        ncomp = optimal.ncomp,
                        keepX = optimal.keepX)

#sPLSDA
png(filename=paste0(work_dir,"/sPLSDA_final_comp1.png"), width=60, height=50, units="cm", res=200)
print({plotIndiv(y_uniq.final.splsda, comp = c(1,1),
        group = sample.group, ind.names = FALSE,
       ellipse = TRUE,
       legend = TRUE,
       title=paste0("sPLS-DA, comp ","1"),
       size.title=rel(3),
       size.xlabel=rel(2), size.ylabel=rel(2),
       size.axis=rel(0.8),
       size.legend.title=rel(2.3), size.legend=rel(2)
       )})
dev.off()


#Heatmap
legend=list(legend=c("Likely NDR","Likely DR"), col=c("#F68B33","#4AA8D8"), title="Likely dr & ndr", cex=1.4)
png(filename=paste0(work_dir,"sPLSDA_final_heatmap_comp1.png"), width=70, height=30, units="cm", res=200)
sample.col <- sample.group
sample.col <- gsub("Likely NDR","#F68B33",sample.col)
sample.col <- gsub("Likely DR","#4AA8D8",sample.col)
print({cim(y_uniq.final.splsda, row.sideColors = sample.col,
            keysize=c(0.2,0.9),
            legend = legend)
    })
dev.off()

#Loading Plot
png(filename=paste0(work_dir,"sPLSDA_loading_Plot.png"),width=30, height=30, units="cm", res=200)
print({
plotLoadings(y_uniq.final.splsda, comp=1, method = 'median', contrib = 'max', size.name=0.7, ndisplay=20, title="Loading Plot comp1")
})
dev.off()

