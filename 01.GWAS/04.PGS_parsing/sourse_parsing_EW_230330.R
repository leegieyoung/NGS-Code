#EW
EWparsing <- function(Start, End){
for (i in (Start:End)){
print(paste0("Num : ", i)) #Check number
#test <- read.table(gzfile(paste0(dir,"PGS000001.txt.gz")), sep='\t')
test <- read.table(gzfile(paste0(dir,EWlist[i],".txt.gz")), sep='\t')
rsID <- grep("rsID",test[1,])
chr <- grep("chr_name",test[1,])
bp <- grep("chr_position",test[1,])
ea <- grep("effect_allele",test[1,])
EW <- grep("effect_weight",test[1,])
ew <- grep("effect_weight", test[1,])

oO <- order(test[,EW], decreasing=T)
test <- test[oO,]

fi <- cbind(test[,rsID], test[,chr], test[,bp], test[,ea], test[,EW], test[,ew])
colnames(fi) <- fi[1,]
fi <- as.data.frame(fi)
fi <- fi[-1,];

#ASA chip and Impute
wA_EW <- which(fi[,rsID] %in% ASAsnplist)
wI_EW <- which(fi[,rsID] %in% IMPUTEsnplist) 
fi_ASA <- as.data.frame(fi[wA_EW,])

#DOI
con <- gzfile(paste0(dir,EWlist[i],".txt.gz"))
if(nrow(na.omit(fi_ASA)) != 0){
fi_ASA$DOI <- grep("^#citation", readLines(con), value=TRUE)}
close(con)
fi_IMP <- fi[wI_EW,]

write.table(fi_ASA, paste0(result_dir,"/EW/",EWlist[i],".A_EW.txt"), sep='\t', quote=F,col.names=T, row.names=F)
write.table(fi_IMP, paste0(result_dir,"/EW/",EWlist[i],".I_EW.txt"), sep='\t', quote=F,col.names=T, row.names=F)

}
}

