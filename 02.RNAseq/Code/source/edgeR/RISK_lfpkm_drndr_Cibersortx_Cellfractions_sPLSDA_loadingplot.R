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

sample.meta <- read.csv("/data/keeyoung/scRNA/iCD/01.Samplelist/RISK_cohort.csv", sep=',')
#Extract
dataset <- grep("RISK",sample.meta$dataset)
sample.meta <- sample.meta[dataset,]
NR  <- grep("^NR$", sample.meta$response)
R <- grep("^R$", sample.meta$response)
sample.meta <- sample.meta[c(NR,R),]

sample.ID <- sample.meta$X
sample.ID <- gsub("-",".", sample.ID)
head(sample.meta)
sample.pheno <- sample.meta$response
names(sample.pheno) <- sample.ID
sample.group <- factor(sample.pheno)

y_uniq <- read.table("/data/keeyoung/scRNA/cibersortx/output/02.Imputed_Cell_Fraction/CD_infla/lfpkm/CIBERSORTx_Adjusted.txt", sep='\t')

y_uniq <- y_uniq[na.omit(match(sample.ID,rownames(y_uniq))),]
dim(y_uniq)

sample.group <- sample.group[na.omit(match(rownames(y_uniq), names(sample.group)))]

print({"Sample group and Y_uniq are same?"})
all.equal(names(sample.group), rownames(y_uniq))
#lfpkm
y_uniq <- y_uniq[,c(-9,-13,-14,-15,-23,-24,-25)]

work_dir <- '/data/keeyoung/gy_RNA/06.output/RISKcohort/04.sPLSDA/lfpkm/'
write.table(y_uniq, paste0(work_dir, "y_uniq.txt"), sep='\t', quote=F)
#Feature select, because of Down sampling
y_uniq <- y_uniq[,c(-13,-18)]

y_uniq_lcpm.splsda <- splsda(y_uniq, sample.group, ncomp=10)

#PLSDA
png(filename=paste0(work_dir,"sPLSDA_plsda.png"), width=60, height=50, units="cm", res=200)
print({plotIndiv(y_uniq_lcpm.splsda, comp = 1:2,
        group = sample.group, ind.names = FALSE,
        ellipse = TRUE,
        legend = TRUE,
        title="PLS-DA"
)})
dev.off()

#Perf
y_uniq_lcpm.perf.splsda <- perf(y_uniq_lcpm.splsda, validation = "Mfold",
                          folds = 5, nrepeat = 100, # use repeated cross-validation
                          progressBar = TRUE, auc = TRUE) # include AUC values
saveRDS(y_uniq_lcpm.perf.splsda, paste0(work_dir,"sPLSDA_perf.rds"))

png(filename=paste0(work_dir,"sPLSDA_perf.png"), width=60, height=50, units="cm", res=200)
print({plot(y_uniq_lcpm.perf.splsda, col = color.mixo(5:7), sd = TRUE,
     legend.position = "horizontal")})
dev.off()

#Check ncomp
y_uniq_lcpm.perf.splsda$choice.ncomp
#ncomp=1

list.keepX <- c(1:10)
y_uniq_lcpm.tune.splsda <- tune.splsda(y_uniq, sample.group, ncomp=1,
                        validation = 'Mfold',
                        folds= 5, nrepeat = 50,
                       dist = 'centroids.dist',
                       measure = 'BER',
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
       title=paste0("sPLS-DA, comp ","& 1,2"),
       size.title=rel(3),
       size.xlabel=rel(2), size.ylabel=rel(2),
       size.axis=rel(0.8),
       size.legend.title=rel(2.3), size.legend=rel(2)
       )})
dev.off()



#Heatmap
legend=list(legend=c("R","NR"), col=c("#F68B33","#C2C2C2"), title="DR & NDR", cex=1.4)
png(filename=paste0(work_dir,"sPLSDA_final_heatmap_comp1.png"), width=70, height=30, units="cm", res=200)
sample.col <- sample.group
sample.col <- gsub("^R$","#F68B33",sample.col)
sample.col <- gsub("^NR$","#C2C2C2",sample.col)
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

