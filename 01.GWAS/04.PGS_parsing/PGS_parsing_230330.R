dir <- "/mnt/nas/PGScatalogScores/"
result_dir <- "/mnt/nas/gylee/0.GWAS/9.etc/PGS/1.Output/"
Input_dir <- "/mnt/nas/gylee/0.GWAS/9.etc/PGS/0.Input/PGS/"
#---
TotalPGSlist <- scan(paste0(Input_dir,"PGS.list"), what=character(0))
ORlist <- scan(paste0(Input_dir,"166OR.list"), what=character(0))
ASAsnplist <- scan(paste0("/mnt/nas/gylee/0.GWAS/9.etc/PGS/SNPlist/","ASAsnp.list"), what=character(0))
IMPUTEsnplist <- scan(paste0("/mnt/nas/gylee/0.GWAS/9.etc/PGS/SNPlist/","IMPUTE.list"), what=character(0))
#OR
op_ASA_list <- c()
op_IMP_list <- c()
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
#top4 & ASAchip 230321-------
#top4 <- head(which(fi[,rsID] %in% ASAsnplist),4)
#fi_top4 <- fi[top4,]

#write.table(fi_top4, paste0(result_dir,"OR/",ORlist[i],".OR.txt"), sep='\t', quote=F,col.names=T)}
#230321-------

#ASAchip 230321--------
wA_OR <- which(fi[,rsID] %in% ASAsnplist)
wI_OR <- which(fi[,rsID] %in% IMPUTEsnplist) 
fi_ASA <- fi[wA_OR,]

#
con <- gzfile(paste0(dir,ORlist[i],".txt.gz"), "rt")
fi_ASA$DOI <- grep("^#citation", readLines(con), value=TRUE)
fi_IMP <- fi[wI_OR,]
write.table(fi_ASA, paste0(result_dir,"/OR/",ORlist[i],".A_OR.txt"), sep='\t', quote=F,col.names=T)
write.table(fi_IMP, paste0(result_dir,"/OR/",ORlist[i],".I_OR.txt"), sep='\t', quote=F,col.names=T)

op_ASA_list <- c(op_ASA_list, fi_ASA[,rsID])
op_IMP_list <- c(op_IMP_list, fi_IMP[,rsID])
}

#Effect_weight
EWlist <- scan(paste0(Input_dir,"1852EW.list"), what=character(0))
for (i in (1:length(EWlist))){
test <- read.table(gzfile(paste0(dir,EWlist[i],".txt.gz")), sep='\t')
rsID <- grep("rsID",test[1,])
chr <- grep("chr_name",test[1,])
bp <- grep("chr_position",test[1,])
ea <- grep("effect_allele",test[1,])
OR <- grep("OR",test[1,])
ew <- grep("effect_weight", test[1,])

fi <- cbind(test[,rsID], test[,chr], test[,bp], test[,ea], test[,OR], test[,ew])
colnames(fi) <- fi[1,]
fi <- fi[-1,];
fiew <- grep("effect_weight", colnames(fi))
fi[,fiew] <- as.numeric(fi[,fiew])
fiewO <- order(fi[,fiew], decreasing=T)
fi <- fi[fiewO,]

#top4 & ASAchip 230321--------
#top4 <- head(which(fi[,rsID] %in% ASAsnplist),4)
#fi_top4 <- fi[top4,]

#write.table(fi_top4, paste0(result_dir,"EW/",EWlist[i],".EW.txt"), sep='\t', quote=F,col.names=T)}
#230321--------
#ASAchip 230321--------
wA_EW <- which(fi[,rsID] %in% ASAsnplist)
wI_EW <- which(fi[,rsID] %in% IMPUTEsnplist)
fi_ASA <- fi[wA_EW,]
fi_IMP <- fi[wI_EW,]
write.table(fi_ASA, paste0(result_dir,"/EW/",EWlist[i],".A_EW.txt"), sep='\t', quote=F,col.names=T)
write.table(fi_IMP, paste0(result_dir,"/EW/",EWlist[i],".I_EW.txt"), sep='\t', quote=F,col.names=T)

op_ASA_list <- c(op_ASA_list, fi_ASA[,rsID])
op_IMP_list <- c(op_IMP_list, fi_IMP[,rsID])
}

write.table(warnings, paste0(result_dir,"EW/","warning_",".EW.txt"), sep='\t', quote=F,col.names=T)


