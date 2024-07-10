#!/bin/Rscript
library(ggplot2);
SAMPLE <- Sys.getenv("SAMPLE")
RESULT_DIR <- Sys.getenv("RESULT_DIR")
print(paste0("Analysis folder : ",RESULT_DIR))
Bhm <- read.table(paste0(RESULT_DIR,"R_AfterQC_check.txt"), header=T, comment.char = "")

Bhm$HET_RATE = (Bhm$OBS_CT - Bhm$O.HOM.)/Bhm$OBS_CT


png(filename=paste0(RESULT_DIR,"R_AfterQC_check.png"), width=20, height=16, units="cm", res=200)
p <- ggplot(Bhm, aes(x=HET_RATE, y=F_MISS, color=Pheno))
p <- p + geom_point()
p <- p + coord_cartesian(xlim = c(min(Bhm$HET_RATE)-0.0001, max(Bhm$HET_RATE)+0.0001), 
	ylim = c(0, max(Bhm$F_MISS)+0.0001))
print({p})
dev.off()


#Ahm <- read.table(paste0(RESULT_DIR,"R_AfterQC_check.txt"), header=T, comment.char = "")
 
#Ahm$HET_RATE = (Ahm$OBS_CT - Ahm$O.HOM.)/Ahm$OBS_CT
 
#print({min(Bhm$F_MISS)})
#print({max(Bhm$F_MISS)})
#print({min(Ahm$F_MISS)})
#print({max(Ahm$F_MISS)})

#png(filename=paste0(RESULT_DIR,"R_AfterQC_check.png"), width=20, height=16, units="cm", res=200)
#p <- ggplot(Ahm, aes(x=HET_RATE, y=F_MISS, color=Pheno))
#p <- p + geom_point()
#p <- p + coord_cartesian(xlim = c(min(Bhm$HET_RATE)-0.0001, max(Bhm$HET_RATE)+0.0001), 
#	ylim = c(0, max(Bhm$F_MISS)+0.0001))
#print({p})
#dev.off()
 

