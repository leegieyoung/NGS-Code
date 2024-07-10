iCD <- read.table("/data/keeyoung/gy_RNA/02.featureCounts/CRC/CRC.txt", sep='\t', head=T)
colnames(iCD) <- gsub(".txt","",colnames(iCD))
pre_dir <- "/data/keeyoung/gy_RNA/06.output/"

sample.meta <- read.table("/data/keeyoung/gy_RNA/01.Samplelist/CRCsamples.txt", sep='\t')
colnames(sample.meta) <- sample.meta[1,]
sample.meta <- sample.meta[-1,]
rownames(sample.meta) <- sample.meta[,1]
sample.meta <- sample.meta[,-1]

sample.ID <- rownames(sample.meta)
sample.pheno <- sample.meta$Class
names(sample.pheno) <- sample.ID
m <- match(names(sample.pheno),colnames(iCD))
iCD <- iCD[,m]

sample.CMS2 <- as.numeric(sample.meta$CMS2)
names(sample.CMS2) <- sample.ID
sample.CMS3 <- as.numeric(sample.meta$CMS3)
names(sample.CMS3) <- sample.ID
sample.CMS4 <- as.numeric(sample.meta$CMS4)
names(sample.CMS4) <- sample.ID
sample.CMS1 <- as.numeric(sample.meta$CMS1)
names(sample.CMS1) <- sample.ID

sample.MSI <- sample.meta$MSI
names(sample.MSI) <- sample.ID

sample.Type <- sample.meta$Type
names(sample.Type) <- sample.ID

#source(paste0("/data/keeyoung/gy_RNA/Code/source/","source_RNAseq_normalized_ensemblToGenesymbol.R")) #Run day is 220919
RNAseq("CRC",iCD)

