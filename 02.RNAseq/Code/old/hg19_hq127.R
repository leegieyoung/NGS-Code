library(tximport)
library(edgeR)
library(mixOmics)
library(org.Hs.eg.db)

#tximport
tx2gene <- read.delim("/data/keeyoung/REFERENCE/RNA/gtf/tx2gene.GRCh37.75.tab")
sample.ID <- scan("/data/keeyoung/gy_RNA/01.Sample/127.IBD.sample", what=character(0))
sample.pheno <- scan("/data/keeyoung/gy_RNA/01.Sample/127.IBD.pheno", what=character(0))
pca.pheno <- scan("/data/keeyoung/gy_RNA/01.Sample/127.IBD.pheno", what=character(0))
names(sample.pheno) <- sample.ID
salmon_root <- "/data/keeyoung/gy_RNA/02.hg19_salmon_result/"
files <- file.path(salmon_root, sample.ID, "quant.sf")
names(files) <- sample.ID

sample.tx <- list()
sample.tx[["salmon"]] <- tximport(files, type="salmon", tx2gene=tx2gene, txOut=FALSE, countsFromAbundance="lengthScaledTPM")

cnt <- sample.tx$salmon$counts

#condition
sample.group <- factor(sample.pheno)
sample.group <- factor(substr(sample.pheno,1,5))
sample.design <- model.matrix(~0 + sample.group)
colnames(sample.design) <- levels(sample.group)

#edgeR
sample.df.cnt <- as.data.frame(sample.tx$salmon$counts)
sample.DGElist <- DGEList(sample.df.cnt, group=sample.group, genes=sample.df.cnt[,0,drop=FALSE])
sample.filterExpr <- sample.DGElist[filterByExpr(sample.DGElist),]
sample.calNorm <- calcNormFactors(sample.filterExpr)
sample.estimateDisp <- estimateDisp(sample.calNorm,sample.design)

#pca
sample.cpm <- cpm(sample.estimateDisp$counts,log=TRUE)
sample.cpm.trans <- t(sample.cpm)
sample.pca <- pca(sample.cpm.trans, ncomp=10, center=TRUE, scale=FALSE)

png(filename="/data/keeyoung/gy_RNA/03.Plot/hg19_127_PCA.png")
plotIndiv(sample.pca, ncomp=c(1,2), ind.names=TRUE, group=pca.pheno, legend=TRUE, guide = "none")
dev.off()

#DEG
DEG_fit <- glmQLFit(sample.estimateDisp, sample.design)
#[1] "CD_in" "CD_no" "UC_in" "UC_no"
CDvsUC_in_DEG <- glmQLFTest(DEG_fit, contrast = c(1,0,-1,0))
CDvsUC_no_DEG <- glmQLFTest(DEG_fit, contrast = c(0,1,0,-1))
INvsNO_CD_DEG <- glmQLFTest(DEG_fit, contrast = c(1,-1,0,0))
INvsNO_UC_DEG <- glmQLFTest(DEG_fit, contrast = c(0,0,1,-1))

CDvsUC_in <- decideTestsDGE(CDvsUC_in_DEG)
CDvsUC_no <- decideTestsDGE(CDvsUC_no_DEG)
INvsNO_CD <- decideTestsDGE(INvsNO_CD_DEG)
INvsNO_UC <- decideTestsDGE(INvsNO_UC_DEG)

CDvsUC_in_FDR <- p.adjust(CDvsUC_in_DEG$table$PValue, method="BH")
CDvsUC_no_FDR <- p.adjust(CDvsUC_no_DEG$table$PValue, method="BH")
INvsNO_CD_FDR <- p.adjust(INvsNO_CD_DEG$table$PValue, method="BH")
INvsNO_UC_FDR <- p.adjust(INvsNO_UC_DEG$table$PValue, method="BH")

CDvsUC_in_result <- cbind(CDvsUC_in_DEG$table, CDvsUC_in_FDR)
CDvsUC_no_result <- cbind(CDvsUC_no_DEG$table, CDvsUC_no_FDR)
INvsNO_CD_result <- cbind(INvsNO_CD_DEG$table, INvsNO_CD_FDR)
INvsNO_UC_result <- cbind(INvsNO_UC_DEG$table, INvsNO_UC_FDR)


write.csv(CDvsUC_in_result, "/data/keeyoung/gy_RNA/04.CSV/hg19_127_CDvsUC_in.csv", row.names=TRUE)
write.csv(CDvsUC_no_result, "/data/keeyoung/gy_RNA/04.CSV/hg19_127_CDvsUC_no.csv", row.names=TRUE)
write.csv(INvsNO_CD_result, "/data/keeyoung/gy_RNA/04.CSV/hg19_127_INvsNO_CD.csv", row.names=TRUE)
write.csv(INvsNO_UC_result, "/data/keeyoung/gy_RNA/04.CSV/127_INvsNO_UC.csv", row.names=TRUE)





