data<- read.table(file="merge_control_case_NoNA_genome_MDS.mds",header=TRUE)
race<- read.table(file="phenofile.txt",header=TRUE)
datafile<- merge(data,race,by=c("IID","FID"))
head(datafile)

pdf("MDS.pdf",width=7,height=7)
for (i in 1:nrow(datafile))
{
if (datafile[i,14]=="Control") {plot(datafile[i,4],datafile[i,5],type="p",xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),xlab="MDS Component 1",ylab="MDS Component 2",pch=1,cex=0.5,col="green")}
par(new=T)
if (datafile[i,14]=="Case") {plot(datafile[i,4],datafile[i,5],type="p",xlim=c(-0.1,0.2),ylim=c(-0.15,0.1),xlab="MDS Component 1",ylab="MDS Component 2",pch=3,cex=0.7,col="black")}
par(new=T)
}

abline(v=-0.035,lty=3)
abline(h=0.035,lty=3)
legend("topright", pch=c(1,3),c("Control","Case"),col=c("green","black"),bty="o",cex=1)
