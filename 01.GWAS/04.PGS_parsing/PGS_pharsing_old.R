dir <- "/mnt/data/PGScatalogScores/"
result_dir <- "/mnt/data/PGScatalogScores/1.Output/"
#---
TotalPGSlist <- scan("/mnt/data/PGScatalogScores/PGS.list", what=character(0))
ORlist <- scan("/mnt/data/PGScatalogScores/1.Output/rsID.OR.list/166OR.list", what=character(0))
ASAsnplist <- scan("/mnt/data/gylee/ASAsnplist/ASAsnpList.txt", what=character(0))

#OR
for (i in (1:length(ORlist))){
#test <- read.table(gzfile(paste0(dir,"PGS000001.txt.gz")), sep='\t')
test <- read.table(gzfile(paste0(dir,ORlist[i],".txt.gz")), sep='\t')
rsID <- grep("rsID",test[1,])
chr <- grep("chr_name",test[1,])
bp <- grep("chr_position",test[1,])
ea <- grep("effect_allele",test[1,])
OR <- grep("OR",test[1,])
ew <- grep("effect_weight", test[1,])

oO <- order(test[,OR], decreasing=T)
test <- test[oO,]

fi <- cbind(test[,rsID], test[,chr], test[,bp], test[,ea], test[,OR], test[,ew])
colnames(fi) <- fi[1,]
fi <- fi[-1,];

write.table(fi, paste0(result_dir,"OR/",ORlist[i],".OR.txt"), sep='\t', quote=F,col.names=T)
}

#Effect_weight
EWlist <- scan("/mnt/data/PGScatalogScores/1.Output/rsID.ea.list/1852EW.list", what=character(0))
for (i in (1:length(EWlist))){
test <- read.table(gzfile(paste0(dir,EWlist[i],".txt.gz")), sep='\t')
rsID <- grep("rsID",test[1,])
chr <- grep("chr_name",test[1,])
bp <- grep("chr_position",test[1,])
ea <- grep("effect_allele",test[1,])
OR <- grep("OR",test[1,])
ew <- grep("effect_weight", test[1,])

ewO <- order(test[,ew], decreasing=T)
test <- test[ewO,]
fi <- cbind(test[,rsID], test[,chr], test[,bp], test[,ea], test[,OR], test[,ew])
colnames(fi) <- fi[1,]
fi <- fi[-1,];
write.table(fi, paste0(result_dir,"EW/",EWlist[i],".EW.txt"), sep='\t', quote=F,col.names=T)
}

