library(tximport)
library(edgeR)
library(mixOmics)
library(org.Hs.eg.db)

#
sample.ID <- scan("/data/keeyoung/gy_RNA/01.Sample/127.IBD.sample", what=character(0))
sample.pheno <- scan("/data/keeyoung/gy_RNA/01.Sample/127.IBD.pheno", what=character(0))
pca.pheno <- scan("/data/keeyoung/gy_RNA/01.Sample/127.IBD.pheno", what=character(0))
names(sample.pheno) <- sample.ID

#condition - record - Pheno
sample.group <- factor(sample.pheno)
sample.group <- factor(substr(sample.pheno,1,5))
sample.design <- model.matrix(~0 + sample.group)
colnames(sample.design) <- levels(sample.group)

#hisat-featurecounts result load
fc127 <- read.csv("/data/keeyoung/gy_RNA/02.featureCounts/gene_name/127_P_M_comma.scv", sep=',', head=T)
fc127 <- fc127[,-1]
fc127 <- t(fc127)

#pheno : dot -> hyphen
fc127.names <- scan("/data/keeyoung/gy_RNA/02.featureCounts/gene_name/127.IBD.hyphen.list", what=character(0))
fc127.ensembl_id <- scan("/data/keeyoung/REFERENCE/RNA/gtf/60650.txt",what=character(0))
colnames(fc127) <- fc127.names
rownames(fc127) <- fc127.ensembl_id

#edgeR
cnt <- fc127
sample.DGElist <- DGEList(cnt, group=sample.group, genes=cnt[,0,drop=FALSE])
sample.filterExpr <- sample.DGElist[filterByExpr(sample.DGElist),]
sample.calNorm <- calcNormFactors(sample.filterExpr)
sample.estimateDisp <- estimateDisp(sample.calNorm,sample.design)

#pca
sample.cpm <- cpm(sample.estimateDisp$counts,log=TRUE)
sample.cpm.trans <- t(sample.cpm)
sample.pca <- pca(sample.cpm.trans, ncomp=10, center=TRUE, scale=FALSE)

png(filename="/data/keeyoung/gy_RNA/03.Plot/fc127_PCA.png")
plotIndiv(sample.pca, ncomp=c(1,2), ind.names=TRUE, group=pca.pheno, legend=TRUE, guide = "none")
dev.off()

#ensembl -> ENTREZID
sample.DGElist$genes$Symbol <- mapIds(org.Hs.eg.db, keys=row.names(sample.DGElist), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
y <- sample.DGElist
y <- y[filterByExpr(y, group = sample.group),]
y <- calcNormFactors(y)
y <- estimateDisp(y,sample.design, robust=TRUE)

y.cpm <- cpm(y$counts,log=TRUE)
y.cpm.trans <- t(y.cpm)
y.pca <- pca(y.cpm.trans, ncomp=10, center=TRUE, scale=FALSE)

#DEG
#[1] "CD_in" "CD_no" "UC_in" "UC_no"
DEG_fit <- glmQLFit(y, sample.design)

CDvsUC_in_DEG <- glmQLFTest(DEG_fit, contrast = c(1,0,-1,0))
CDvsUC_no_DEG <- glmQLFTest(DEG_fit, contrast = c(0,1,0,-1))
INvsNO_CD_DEG <- glmQLFTest(DEG_fit, contrast = c(1,-1,0,0))
INvsNO_UC_DEG <- glmQLFTest(DEG_fit, contrast = c(0,0,1,-1))

CDvsUC_in <- decideTestsDGE(CDvsUC_in_DEG)
CDvsUC_no <- decideTestsDGE(CDvsUC_no_DEG)
INvsNO_CD <- decideTestsDGE(INvsNO_CD_DEG)
INvsNO_UC <- decideTestsDGE(INvsNO_UC_DEG)

##Calc-FDR
CDvsUC_in_FDR <- p.adjust(CDvsUC_in_DEG$table$PValue, method="BH")
CDvsUC_no_FDR <- p.adjust(CDvsUC_no_DEG$table$PValue, method="BH")
INvsNO_CD_FDR <- p.adjust(INvsNO_CD_DEG$table$PValue, method="BH")
INvsNO_UC_FDR <- p.adjust(INvsNO_UC_DEG$table$PValue, method="BH")

CDvsUC_in_result <- cbind(CDvsUC_in_DEG$table, CDvsUC_in_FDR)
CDvsUC_no_result <- cbind(CDvsUC_no_DEG$table, CDvsUC_no_FDR)
INvsNO_CD_result <- cbind(INvsNO_CD_DEG$table, INvsNO_CD_FDR)
INvsNO_UC_result <- cbind(INvsNO_UC_DEG$table, INvsNO_UC_FDR)


write.table(CDvsUC_in_result, "/data/keeyoung/gy_RNA/04.CSV/fc127_P_M_CDvsUC_in_ensembl.csv", row.names=TRUE, quote=F, sep=",")
write.table(CDvsUC_no_result, "/data/keeyoung/gy_RNA/04.CSV/fc127_P_M_CDvsUC_no_ensembl.csv", row.names=TRUE, quote=F, sep=",")
write.table(INvsNO_CD_result, "/data/keeyoung/gy_RNA/04.CSV/fc127_P_M_INvsNO_CD_ensembl.csv", row.names=TRUE, quote=F, sep=",")
write.table(INvsNO_UC_result, "/data/keeyoung/gy_RNA/04.CSV/fc127_P_M_INvsNO_UC_ensembl.csv", row.names=TRUE, quote=F, sep=",")

ensembl_cnt <- y$counts
rownames(ensembl_cnt) <- y$genes$Symbol
entrez_cnt <- subset(ensembl_cnt, rownames(ensembl_cnt)!="NA")
entrez_cnt <- DGEList(entrez_cnt, group=sample.group, genes=entrez_cnt[,0,drop=FALSE])
entrez_cnt <- entrez_cnt[filterByExpr(entrez_cnt),]
entrez_cnt <- calcNormFactors(entrez_cnt)
entrez_cnt <- estimateDisp(entrez_cnt,sample.design)

#DEG
#[1] "CD_in" "CD_no" "UC_in" "UC_no"
DEG_fit <- glmQLFit(entrez_cnt, sample.design)

CDvsUC_in_DEG <- glmQLFTest(DEG_fit, contrast = c(1,0,-1,0))
CDvsUC_no_DEG <- glmQLFTest(DEG_fit, contrast = c(0,1,0,-1))
INvsNO_CD_DEG <- glmQLFTest(DEG_fit, contrast = c(1,-1,0,0))
INvsNO_UC_DEG <- glmQLFTest(DEG_fit, contrast = c(0,0,1,-1))

CDvsUC_in <- decideTestsDGE(CDvsUC_in_DEG)
CDvsUC_no <- decideTestsDGE(CDvsUC_no_DEG)
INvsNO_CD <- decideTestsDGE(INvsNO_CD_DEG)
INvsNO_UC <- decideTestsDGE(INvsNO_UC_DEG)

##Calc-FDR
CDvsUC_in_FDR <- p.adjust(CDvsUC_in_DEG$table$PValue, method="BH")
CDvsUC_no_FDR <- p.adjust(CDvsUC_no_DEG$table$PValue, method="BH")
INvsNO_CD_FDR <- p.adjust(INvsNO_CD_DEG$table$PValue, method="BH")
INvsNO_UC_FDR <- p.adjust(INvsNO_UC_DEG$table$PValue, method="BH")

CDvsUC_in_result <- cbind(CDvsUC_in_DEG$table, CDvsUC_in_FDR)
CDvsUC_no_result <- cbind(CDvsUC_no_DEG$table, CDvsUC_no_FDR)
INvsNO_CD_result <- cbind(INvsNO_CD_DEG$table, INvsNO_CD_FDR)
INvsNO_UC_result <- cbind(INvsNO_UC_DEG$table, INvsNO_UC_FDR)


write.table(CDvsUC_in_result, "/data/keeyoung/gy_RNA/04.CSV/fc127_P_M_CDvsUC_in_ENTREZID.csv", row.names=TRUE, quote=F, sep=",")
write.table(CDvsUC_no_result, "/data/keeyoung/gy_RNA/04.CSV/fc127_P_M_CDvsUC_no_ENTREZID.csv", row.names=TRUE, quote=F, sep=",")
write.table(INvsNO_CD_result, "/data/keeyoung/gy_RNA/04.CSV/fc127_P_M_INvsNO_CD_ENTREZID.csv", row.names=TRUE, quote=F, sep=",")
write.table(INvsNO_UC_result, "/data/keeyoung/gy_RNA/04.CSV/fc127_P_M_INvsNO_UC_ENTREZID.csv", row.names=TRUE, quote=F, sep=",")

