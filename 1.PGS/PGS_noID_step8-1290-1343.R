dir <- "/mnt/nas/PGScatalogScores/"
result_dir <- "/mnt/nas/gylee/0.GWAS/9.etc/PGS/1.Output/"
Input_dir <- "/mnt/nas/gylee/0.GWAS/9.etc/PGS/0.Input/PGS/"
#---
#TotalPGSlist <- scan(paste0(Input_dir,"PGS.list"), what=character(0))
#ASAsnplist <- scan(paste0("/mnt/nas/gylee/0.GWAS/9.etc/PGS/SNPlist/","ASAsnp.list"), what=character(0))
noIDlist <- scan(paste0(Input_dir,"1330GRChA1A2.list"), what=character(0))
#IMPUTEsnplist <- scan(paste0("/mnt/nas/gylee/0.GWAS/9.etc/PGS/SNPlist/","IMPUTE.list"), what=character(0))

source(paste0("/mnt/nas/gylee/0.GWAS/Code/1.PGS/","sourse_parsing_noID_OR_230405.R"))
noIDparsing(1290, 1343)
