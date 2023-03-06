raw_data <- read.csv("/data/keeyoung/gy_RNA/02.featureCounts/454.IBD.fc_P_M.merge", sep='\t', head=T)
genelist <- raw_data[,1]
raw_data <- raw_data[,-1]
rownames(raw_data) <- genelist
raw_data <- t(raw_data)
write.table(raw_data, file="/data/keeyoung/gy_RNA/02.featureCounts/trans_454.IBD.fc_P_M.merge", quote=F, sep="\t")



