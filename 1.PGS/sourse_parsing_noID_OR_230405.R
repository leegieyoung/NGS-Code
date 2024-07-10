#noID
noIDparsing <- function(Start, End){
for (i in (Start:End)){
print(paste0("Num : ", i)) #Check number
#test <- read.table(gzfile(paste0(dir,"PGS000001.txt.gz")), sep='\t')
test <- read.table(gzfile(paste0(dir,noIDlist[i],".txt.gz")), sep='\t')
chr <- grep("chr_name",test[1,])
bp <- grep("chr_position",test[1,])
ea <- grep("effect_allele",test[1,])
oa <- grep("other_allele",test[1,])
OR <- grep("OR",test[1,])
ew <- grep("effect_weight", test[1,])

oEW <- order(test[,ew], decreasing=T)
test <- test[oEW,]

fi <- cbind(test[,chr], test[,bp], test[,ea], test[,oa], test[,ew], test[,OR])
colnames(fi) <- fi[1,]
fi <- as.data.frame(fi)
fi <- fi[-1,];

#DOI
con <- gzfile(paste0(dir,noIDlist[i],".txt.gz"))
if(nrow(na.omit(fi)) != 0){
fi$DOI <- grep("^#citation", readLines(con), value=TRUE)}
close(con)

write.table(fi, paste0(result_dir,"/noID/",noIDlist[i],".noID.txt"), sep='\t', quote=F,col.names=T, row.names=F)

}
}

