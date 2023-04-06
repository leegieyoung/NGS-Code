dir <- "/mnt/nas/PGScatalogScores/"
result_dir <- "/mnt/nas/gylee/0.GWAS/9.etc/PGS/1.Output/"
Input_dir <- "/mnt/nas/gylee/0.GWAS/9.etc/PGS/0.Input/PGS/"
#---
TotalPGSlist <- scan(paste0(Input_dir,"PGS.list"), what=character(0))
EWlist <- scan(paste0(Input_dir,"1852EW.list"), what=character(0))
ASAsnplist <- scan(paste0("/mnt/nas/gylee/0.GWAS/9.etc/PGS/SNPlist/","ASAsnp.list"), what=character(0))
IMPUTEsnplist <- scan(paste0("/mnt/nas/gylee/0.GWAS/9.etc/PGS/SNPlist/","IMPUTE.list"), what=character(0))

source(paste0("/mnt/nas/gylee/0.GWAS/Code/1.PGS/","sourse_parsing_EW_230330.R"))
EWparsing(465,696)
