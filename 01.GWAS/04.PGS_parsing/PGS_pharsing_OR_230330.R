dir <- "/mnt/nas/PGScatalogScores/"
result_dir <- "/mnt/nas/gylee/0.GWAS/9.etc/PGS/1.Output/"
Input_dir <- "/mnt/nas/gylee/0.GWAS/9.etc/PGS/0.Input/PGS/"
#---
TotalPGSlist <- scan(paste0(Input_dir,"PGS.list"), what=character(0))
ORlist <- scan(paste0(Input_dir,"166OR.list"), what=character(0))
ASAsnplist <- scan(paste0("/mnt/nas/gylee/0.GWAS/9.etc/PGS/SNPlist/","ASAsnp.list"), what=character(0))
IMPUTEsnplist <- scan(paste0("/mnt/nas/gylee/0.GWAS/9.etc/PGS/SNPlist/","IMPUTE.list"), what=character(0))
#OR
for (i in (1:length(ORlist))){
print(paste0("Num : ", i)) #Check number
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
fi <- as.data.frame(fi)

#ASA chip and Impute
wA_OR <- which(fi[,rsID] %in% ASAsnplist)
wI_OR <- which(fi[,rsID] %in% IMPUTEsnplist) 
fi_ASA <- fi[wA_OR,]

#DOI
con <- gzfile(paste0(dir,ORlist[i],".txt.gz"), "rt")
fi_ASA$DOI <- grep("^#citation", readLines(con), value=TRUE)
fi_IMP <- fi[wI_OR,]

write.table(fi_ASA, paste0(result_dir,"/OR/",ORlist[i],".A_OR.txt"), sep='\t', quote=F,col.names=T)
write.table(fi_IMP, paste0(result_dir,"/OR/",ORlist[i],".I_OR.txt"), sep='\t', quote=F,col.names=T)

}


