library(Rsubread)

gtf="/data/keeyoung/gy_RNA/Code/featurecounts/Homo_sapiens.GRCh38.100.gtf.gz"
pre_dir="/data/keeyoung/gy_RNA/Code/featurecounts/"


FeatureCounts <-function(bam){
	bamfile <-paste0(pre_dir,"hisat2.dta.sorted.bam") 
	fct=featureCounts(file=bamfile, isPairedEnd=TRUE, countMultiMappingReads=TRUE, countReadPairs=TRUE, isGTFAnnotationFile=TRUE, nthreads=30, annot.ext=gtf, tmpDir="/data/keeyoung/gy_RNA/Code/featurecounts/temp")
write.table(fct$counts, file=paste0("/data/keeyoung/gy_RNA/Code/featurecounts/",bam,".txt"), row.names=T, col.names=T, sep=' ', quote=F)
}
