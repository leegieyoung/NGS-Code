library(Seurat)
read10x <- function(a,dir) {
    names <- c("barcodes.tsv", "genes.tsv", "matrix.mtx")
    for (i in 1:length(names)) {
        R.utils::gunzip(paste0(dir, "/",a,"_" ,names[i], ".gz"))
	pre <- paste0(dir, "/",a ,"_" ,names[i])
	new <- paste0(dir, "/",names[i])
	file.rename(from = pre , to =new)
    }
#   mat <- Seurat::Read10X(dir)
}



